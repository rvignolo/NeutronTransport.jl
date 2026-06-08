# Numeric golden references copied from OpenMOC's committed regression outputs.
#
# Homogeneous infinite medium:
# https://github.com/mit-crpg/OpenMOC/blob/develop/tests/test_forward_hom_inf_medium/results_true.dat
#
# Pin-cell reference:
# https://github.com/mit-crpg/OpenMOC/blob/develop/tests/test_forward_pin_cell/results_true.dat

const OPENMOC_HOMOGENEOUS_REFERENCE = (
    iterations = 108,
    keff = 1.72307,
    flux_group_1 = 76.86901,
    flux_group_2 = 48.87577,
)

const OPENMOC_PIN_CELL_REFERENCE = (
    iterations = 261,
    keff = 1.04666,
    fluxes = [
        0.6428760,
        1.099273,
        0.5710249,
        0.2507108,
        0.1962766,
        0.4993489,
        1.287533,
        1.104038,
        1.411933,
        0.5559559,
        0.2268887,
        0.1900581,
        0.4444841,
        0.9440460,
    ],
)

function openmoc_homogeneous_reference()
    ref = OPENMOC_HOMOGENEOUS_REFERENCE
    return (
        iterations = ref.iterations,
        keff = ref.keff,
        flux_ratio = ref.flux_group_2 / ref.flux_group_1,
    )
end

function openmoc_pin_cell_reference()
    return OPENMOC_PIN_CELL_REFERENCE
end
