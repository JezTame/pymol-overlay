"""overlay.py -- strict least-squares superposition in PyMOL.

Load with: run overlay.py
Usage: overlay mobile, target [, mode=strict|ordered [, report=auto|0|filename]]

Protein: CA pairs; ligand: named heavy-atom pairs. Single or multiple protein chains; one ligand residue.
Separate objects, first coordinate state only. Moves entire mobile object.
"""

from collections import Counter
from datetime import datetime
from pathlib import Path
import re

import numpy as np
from pymol import cmd


class OverlayError(ValueError):
    pass


def _residue_sort_key(resi):
    """Sort integer PDB residue numbers, then optional insertion codes."""
    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", resi.strip())
    if not m:
        raise OverlayError("Unsupported residue identifier {!r}".format(resi))
    return int(m.group(1)), m.group(2)


def _inspect(selection, side):
    expr = "({})".format(selection)
    atoms = cmd.get_model(expr, state=1).atom
    if not atoms:
        raise OverlayError("{} selection contains no atoms: {}".format(side, selection))

    objects = set(cmd.get_object_list(expr))
    if len(objects) != 1:
        raise OverlayError("{} must select atoms from exactly one object (found {}).".format(
            side, sorted(objects)))
    obj = next(iter(objects))
    if cmd.get_type(obj) != 'object:molecule':
        raise OverlayError("{} is not a molecular object.".format(obj))
    if cmd.count_states(obj) != 1:
        raise OverlayError("{} has multiple states; this draft requires one.".format(obj))

    # Define protein using PyMOL's polymer classification, not CA alone.
    protein_atoms = cmd.get_model("({}) and polymer.protein".format(selection), state=1).atom
    if protein_atoms:
        # A protein selection may span several chains of the same parent object.
        # Chain identity is (segment, chain), so chains with repeated names in
        # different segments remain distinct.
        ca = [a for a in protein_atoms if a.name == 'CA']
        if len(ca) < 3:
            raise OverlayError("{} protein has fewer than 3 C-alpha atoms.".format(side))
        keys = [(a.segi, a.chain, a.resi) for a in ca]
        duplicates = [k for k, n in Counter(keys).items() if n != 1]
        if duplicates:
            raise OverlayError("{} has duplicate C-alpha identifiers (possibly alternate locations): {}".format(
                side, duplicates[:5]))
        ca.sort(key=lambda a: (a.segi, a.chain, _residue_sort_key(a.resi)))
        return obj, 'protein', ca

    # Nonprotein: require one chemically defined residue, not water/ions.
    if cmd.count_atoms("({}) and (polymer.nucleic or solvent)".format(selection)):
        raise OverlayError("{} is not a protein or a supported single ligand.".format(side))
    heavy = [a for a in atoms if a.symbol.upper() not in ('H', 'D')]
    if len(heavy) < 3:
        raise OverlayError("{} ligand has fewer than 3 heavy atoms.".format(side))
    residues = {(a.segi, a.chain, a.resi, a.resn) for a in heavy}
    if len(residues) != 1:
        raise OverlayError("{} ligand selection must contain exactly one residue; found {}.".format(
            side, len(residues)))
    duplicates = [name for name, n in Counter(a.name for a in heavy).items() if n != 1]
    if duplicates:
        raise OverlayError("{} has duplicate ligand atom names: {}".format(side, duplicates))
    heavy.sort(key=lambda a: a.name)
    return obj, 'ligand', heavy


def _correspond(moving, fixed, molecule_type, mode):
    if molecule_type == 'protein':
        mchains = {(a.segi, a.chain) for a in moving}
        fchains = {(a.segi, a.chain) for a in fixed}
        if len(mchains) != len(fchains):
            raise OverlayError("Different numbers of selected protein chains: mobile={}, target={}.".format(
                len(mchains), len(fchains)))

        if len(mchains) == 1:
            # One-chain convenience: chain letters may differ.
            if mode == 'strict':
                mov_ids = [a.resi for a in moving]
                fix_ids = [a.resi for a in fixed]
                if mov_ids != fix_ids:
                    missing = sorted(set(fix_ids) - set(mov_ids), key=_residue_sort_key)
                    extra = sorted(set(mov_ids) - set(fix_ids), key=_residue_sort_key)
                    raise OverlayError(
                        "C-alpha residue numbering mismatch. Missing in mobile: {}; "
                        "extra in mobile: {}. Use mode=ordered only if positions correspond.".format(
                            missing[:12], extra[:12]))
            elif mode == 'ordered':
                if len(moving) != len(fixed):
                    raise OverlayError("Ordered C-alpha matching requires equal counts: mobile={}, target={}.".format(
                        len(moving), len(fixed)))
                print("OVERLAY WARNING: ordered mode does not verify sequence equivalence.")
            else:
                raise OverlayError("Unknown mode {!r}; use strict or ordered.".format(mode))
        else:
            # Multichain safety: never pair different chains merely because
            # their residues happen to share an index or atom order.
            if mode != 'strict':
                raise OverlayError("Multichain fitting currently requires mode=strict.")
            if mchains != fchains:
                raise OverlayError(
                    "Multichain fitting requires identical chain and segment identifiers. "
                    "Mobile chains: {}; target chains: {}. "
                    "Select/rename corresponding chains before fitting.".format(
                        sorted(mchains), sorted(fchains)))
            mkeys = {(a.segi, a.chain, a.resi) for a in moving}
            fkeys = {(a.segi, a.chain, a.resi) for a in fixed}
            if mkeys != fkeys:
                missing = sorted(fkeys - mkeys, key=lambda k: (k[0], k[1], _residue_sort_key(k[2])))
                extra = sorted(mkeys - fkeys, key=lambda k: (k[0], k[1], _residue_sort_key(k[2])))
                raise OverlayError(
                    "Multichain C-alpha correspondence mismatch. "
                    "Missing in mobile: {}; extra in mobile: {}.".format(
                        missing[:12], extra[:12]))
    else:
        if mode != 'strict':
            raise OverlayError("Ligands always use strict atom-name matching (mode=strict).")
        m = {a.name: a for a in moving}
        f = {a.name: a for a in fixed}
        if m.keys() != f.keys():
            missing = sorted(f.keys() - m.keys())
            extra = sorted(m.keys() - f.keys())
            raise OverlayError("Ligand names differ. Missing in mobile: {}; extra in mobile: {}".format(
                missing, extra))
        for name in sorted(m):
            if m[name].symbol.upper() != f[name].symbol.upper():
                raise OverlayError("Element mismatch for {}: mobile={}, target={}".format(
                    name, m[name].symbol, f[name].symbol))

    if len(moving) != len(fixed):
        raise OverlayError("Atom-count mismatch: mobile={}, target={}.".format(len(moving), len(fixed)))
    return np.array([a.coord for a in moving], dtype=float), np.array(
        [a.coord for a in fixed], dtype=float)


def _kabsch(mobile_xyz, target_xyz):
    """Return R, t for target ~= R @ mobile + t, plus post-fit RMSD."""
    cm = mobile_xyz.mean(axis=0)
    ct = target_xyz.mean(axis=0)
    h = (mobile_xyz - cm).T @ (target_xyz - ct)
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = ct - r @ cm
    fitted = mobile_xyz @ r.T + t
    rmsd = np.sqrt(np.mean(np.sum((fitted - target_xyz) ** 2, axis=1)))
    return r, t, float(rmsd)



def _screw_parameters(r, t):
    """Return screw rotation angle, unit axis and signed axial displacement.

    Axis is undefined for a zero-angle rotation. In that case, the
    transformation is a pure translation of magnitude ||t||.
    """
    cos_theta = float(np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0))
    theta = float(np.arccos(cos_theta))
    # Pure translation: any direction along t can serve as an axis.
    if theta < 1e-6:
        magnitude = float(np.linalg.norm(t))
        direction = t / magnitude if magnitude > 1e-12 else None
        return 0.0, direction, magnitude, True

    # An eigenvector of R with eigenvalue 1 lies along the screw axis.
    # Works at 180 degrees too, unlike division by sin(theta).
    eigenvalues, eigenvectors = np.linalg.eig(r)
    k = int(np.argmin(np.abs(eigenvalues - 1.0)))
    direction = np.real(eigenvectors[:, k])
    direction /= np.linalg.norm(direction)

    # Orient the eigenvector consistently with the rotation's handedness.
    # At exactly 180 degrees the sign convention is inherently ambiguous.
    skew = np.array([r[2, 1]-r[1, 2],
                     r[0, 2]-r[2, 0],
                     r[1, 0]-r[0, 1]])
    if np.linalg.norm(skew) > 1e-8 and np.dot(direction, skew) < 0:
        direction = -direction
    elif np.linalg.norm(skew) <= 1e-8:
        if direction[int(np.argmax(np.abs(direction)))] < 0:
            direction = -direction
    axial = float(np.dot(direction, t))
    return float(np.degrees(theta)), direction, axial, False


def _report_path(report, mo, fo):
    """Resolve report option; 0/off/none disables reporting."""
    value = str(report).strip()
    if value.lower() in ('0', 'off', 'none', 'no', 'false'):
        return None
    if value.lower() in ('auto', '1', 'yes', 'true'):
        safe = lambda name: re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('._') or 'object'
        return Path('overlay_{}_onto_{}.txt'.format(safe(mo), safe(fo)))
    return Path(value).expanduser()




def _arrow_commands(r, t, moving_xyz, axis, angle, axial):
    """Paste-ready cgo_arrow commands illustrating the fitted rigid motion.

    The reported COM is the unweighted centroid of the matched fitting atoms,
    not the molecular centre of mass. All coordinates are in the original
    target coordinate frame, prior to applying the transformation in PyMOL.
    """
    before = moving_xyz.mean(axis=0)
    after = r @ before + t
    commands = []

    def add(start, end, color, name):
        if np.linalg.norm(end - start) < 1e-4:
            return
        color_str = color
        fmt = lambda v: '[' + ', '.join('{:.4f}'.format(float(x)) for x in v) + ']'
        commands.append('cgo_arrow {}, {}, 0.25, 0, -1, -1, {}, {}'.format(
            fmt(start), fmt(end), color_str, repr(name)))

    add(before, after, 'red', 'COM_displacement')
    if axis is None or angle < 1e-4:
        return commands, 'No unique rotation axis (zero or negligible rotation).'

    # A point on the screw axis satisfies (I-R) p = t - d*axis,
    # where d is signed screw displacement. Choose its position along
    # the axis so it lies near the midpoint of the fitting-atom centroids.
    midpoint = 0.5 * (before + after)
    p0 = np.linalg.lstsq(np.eye(3) - r, t - axial * axis, rcond=None)[0]
    p = p0 + axis * np.dot(midpoint - p0, axis)
    proj_before = p + axis * np.dot(before - p, axis)
    proj_after = p + axis * np.dot(after - p, axis)

    span = max(5.0, float(np.max(np.linalg.norm(moving_xyz - before, axis=1))))
    add(p - axis * span, p + axis * span, 'yellow', 'screw_axis')
    add(before, proj_before, 'cyan', 'COM_before_to_axis')
    add(proj_after, after, 'magenta', 'axis_to_COM_after')
    add(proj_before, proj_after, 'orange', 'axial_shift')
    add(p, midpoint, 'green', 'axis_to_midpoint')
    return commands, None

def _write_report(path, mo, fo, mobile, target, mode, molecule_type, n, rmsd, r, t, moving_xyz):
    angle, axis, axial, pure_translation = _screw_parameters(r, t)
    with path.open('w', encoding='utf-8') as out:
        out.write('PyMOL OVERLAY -- least-squares superposition report\n')
        out.write('Generated: {}\n\n'.format(datetime.now().astimezone().isoformat(timespec='seconds')))
        out.write('Mobile selection: {}\nTarget selection: {}\n'.format(mobile, target))
        out.write('Mobile object: {}\nTarget object: {}\n'.format(mo, fo))
        out.write('Molecule type: {}\nMatching mode: {}\n'.format(molecule_type, mode))
        out.write('Matched atom pairs: {}\nRMSD after fit: {:.6f} Angstrom\n\n'.format(n, rmsd))
        out.write('Coordinate convention: x_target = R @ x_mobile + t\n')
        out.write('Rotation matrix R (3 x 3):\n')
        for row in r:
            out.write('  ' + ' '.join('{: .9f}'.format(float(v)) for v in row) + '\n')
        out.write('Translation vector t (Angstrom):\n  ' +
                  ' '.join('{: .9f}'.format(float(v)) for v in t) + '\n\n')
        out.write('Homogeneous RT matrix (4 x 4, column-vector convention):\n')
        rt = np.eye(4)
        rt[:3, :3] = r
        rt[:3, 3] = t
        for row in rt:
            out.write('  ' + ' '.join('{: .9f}'.format(float(v)) for v in row) + '\n')
        out.write('\nRotation angle: {:.6f} degrees\n'.format(angle))
        if axis is None:
            out.write('Screw axis direction: undefined (identity transformation)\n')
        else:
            out.write('Screw axis direction: ' +
                      ' '.join('{: .9f}'.format(float(v)) for v in axis) + '\n')
        if pure_translation:
            out.write('Pure-translation distance: {:.6f} Angstrom\n'.format(abs(axial)))
            out.write('Screw axial translation: not uniquely defined at zero rotation\n')
        else:
            out.write('Signed screw axial translation: {:.6f} Angstrom\n'.format(axial))
            out.write('Screw axial translation magnitude: {:.6f} Angstrom\n'.format(abs(axial)))
        out.write('\nCGO ARROWS (requires cgo_arrow.py loaded in PyMOL)\n')
        out.write('Run: run cgo_arrow.py (once per session, if not already loaded)\n')
        out.write('Paste the following commands into the PyMOL command line.\n\n')
        commands, note = _arrow_commands(r, t, moving_xyz, axis, angle, axial)
        if note:
            out.write('# ' + note + '\n')
        if commands:
            for command in commands:
                out.write(command + '\n')
        else:
            out.write('# No nonzero displacement to display.\n')


def overlay(mobile, target, mode='strict', report='auto'):
    """Fit one molecular object onto another using validated atom pairs.

    Usage: overlay mobile, target [, mode=strict|ordered [, report=auto|0|filename]]
    Protein: C-alpha atoms (single or several chains). Ligand: named heavy atoms.
    The selection defines fitting atoms; the entire mobile object moves.
    """
    try:
        mode = str(mode).lower().strip()
        mo, mt, ma = _inspect(mobile, 'Mobile')
        fo, ft, fa = _inspect(target, 'Target')
        if mo == fo:
            raise OverlayError("Mobile and target must be different objects.")
        if mt != ft:
            raise OverlayError("Molecule types differ: mobile={}, target={}.".format(mt, ft))
        moving_xyz, fixed_xyz = _correspond(ma, fa, mt, mode)
        if np.linalg.matrix_rank(moving_xyz - moving_xyz.mean(axis=0)) < 2:
            raise OverlayError("Mobile fitting atoms are collinear or degenerate.")
        if np.linalg.matrix_rank(fixed_xyz - fixed_xyz.mean(axis=0)) < 2:
            raise OverlayError("Target fitting atoms are collinear or degenerate.")
        r, t, rmsd = _kabsch(moving_xyz, fixed_xyz)

        # PyMOL accepts a row-major homogeneous 4x4 transformation matrix.
        transform = np.eye(4)
        transform[:3, :3] = r
        transform[:3, 3] = t
        cmd.transform_object(mo, transform.ravel().tolist(), state=1, homogenous=1)

        print("OVERLAY: {} fitting ({}), {} atom pairs".format(mt, mode, len(ma)))
        if mt == "protein":
            print("  Selected protein chains: {}".format(len({(a.segi, a.chain) for a in ma})))
        print("  Mobile object: {}  ->  Target object: {}".format(mo, fo))
        print("  RMSD after fitting: {:.4f} A".format(rmsd))
        destination = _report_path(report, mo, fo)
        if destination is not None:
            try:
                _write_report(destination, mo, fo, mobile, target, mode,
                              mt, len(ma), rmsd, r, t, moving_xyz)
                print("  Report written: {}".format(destination.resolve()))
            except OSError as exc:
                print("OVERLAY WARNING: fit succeeded, but report could not be written: {}".format(exc))
        return rmsd
    except OverlayError as exc:
        print("OVERLAY ERROR: {}".format(exc))
        return None


cmd.extend('overlay', overlay)
