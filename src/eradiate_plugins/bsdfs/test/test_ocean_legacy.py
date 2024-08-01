import pytest
import drjit as dr
import mitsuba as mi


@pytest.mark.slow
def test01_chi2(variants_vec_backends_once_rgb):
    params = {
        'type': 'oceanic_legacy',
        'channel': 3,
        'visual_type': 0,
        'wavelength': 2.2,
        'wind_speed': 2,
        'wind_direction': 0,
        'chlorinity': 19,
        'pigmentation': 0.3,
        'shininess': 10.0,
    }
    sample_func, pdf_func = mi.chi2.BSDFAdapter("oceanic_legacy", params)

    chi2 = mi.chi2.ChiSquareTest(
        domain=mi.chi2.SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
        ires=16,
        res=201
    )

    assert chi2.run()

    chi2._dump_tables()