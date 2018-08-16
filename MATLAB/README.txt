---Main scripts---
polar_interpolation.m - Interpolates rings into polar coordinates

ringFit2D_main.m - Fit ring image in polar coordinates

---Other scripts---

generate_synth_sample.m - Synthesize HEXD data that vary spatially on a sample

miller_indices_rings.m - Compares experimental setup geometry to theoretical values


---Primary functions---
extract_ring.m
FISTA_Circulant.m
computeAWMV.m
computeHVS.m
computeFitError.m

---Secondary functions---
unshifted_basis_matrix_ft_stack_norm2.m
unshifted_basis_matrix_stack_norm2.m
unshifted_basis_matrix_ft_stack_norm.m
unshifted_basis_matrix_stack_norm.m
unshifted_basis_matrix_ft_stack.m
unshifted_basis_matrix_stack.m
gaussian_basis_wrap_2D_norm2.m
shift2D.m
AtR_ft_2D.m
Ax_ft_2D.m
forceMaskToZero.m
forceMaskTOZeroArray.m
zeroPad.m
onePad.m
load_neighbors_awmv.m
azimuthal_projection.m

---Cluster codes---
Files starting in slurm_, wrap_


Anything else is either old, data specific, or just not very useful