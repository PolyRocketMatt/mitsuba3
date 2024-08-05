#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <drjit/dynamic.h>

#include <typeinfo> // For typeid

NAMESPACE_BEGIN(mitsuba)

// Header content - Potentially move somewhere else later
#ifndef OCEAN_PROPS
#define OCEAN_PROPS

template<typename Float, typename Spectrum>
class OceanProperties {
public:
    MI_IMPORT_TYPES()

    using AzimuthStorage = dr::Array<ScalarFloat, 48>;
    using ZenithStorage = dr::Array<ScalarFloat, 24>;
    using Value = dr::Array<ScalarFloat, 1>;

    OceanProperties() {
        // Effective reflectance of whitecaps (Whitlock et al. 1982)
        std::vector<ScalarFloat> wc_wavelengths = {     0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
                                                        2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1,
                                                        3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0 };
        std::vector<ScalarFloat> wc_data = {    0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210, 0.200, 0.190,
                                                0.175, 0.155, 0.130, 0.080, 0.100, 0.105, 0.100, 0.080, 0.045, 0.055,
                                                0.065, 0.060, 0.055, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };

        // Complex index of refraction of water (Hale & Querry 1973)
        std::vector<ScalarFloat> ior_wavelengths = {    0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.345 ,0.375, 0.400, 0.425 ,
                                                        0.445, 0.475, 0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675,
                                                        0.700, 0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925,
                                                        0.950, 0.975, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 
                                                        2.600, 2.650, 2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, 3.050,
                                                        3.100, 3.150, 3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 3.500, 3.600,
                                                        3.700, 3.800, 3.900, 4.000 };
        std::vector<ScalarFloat> ior_real_data = {  1.369, 1.373, 1.362, 1.354, 1.349, 1.346, 1.343, 1.341, 1.339, 1.338,
                                                    1.337, 1.336, 1.335, 1.334, 1.333, 1.333, 1.332, 1.332, 1.331, 1.331,
                                                    1.331, 1.330, 1.330, 1.330, 1.329, 1.329, 1.329, 1.328, 1.328, 1.328, 
                                                    1.327, 1.327, 1.327, 1.324, 1.321, 1.317, 1.312, 1.306, 1.296, 1.279,
                                                    1.242, 1.219, 1.188, 1.157, 1.142, 1.149, 1.201, 1.292, 1.371, 1.426,
                                                    1.467, 1.483, 1.478, 1.467, 1.450, 1.432, 1.420, 1.410, 1.400, 1.385,
                                                    1.374, 1.364, 1.357, 1.351 };
        std::vector<ScalarFloat> ior_cplx_data = {  1.10e-07, 4.90e-08, 3.35e-08, 2.35e-08, 1.60e-08, 
                                                    1.08e-08, 6.50e-09, 3.50e-09, 1.86e-09, 1.30e-09, 
                                                    1.02e-09, 9.35e-10, 1.00e-09, 1.32e-09, 1.96e-09, 
                                                    3.60e-09, 1.09e-08, 1.39e-08, 1.64e-08, 2.23e-08, 
                                                    3.35e-08, 9.15e-08, 1.56e-07, 1.48e-07, 1.25e-07, 
                                                    1.82e-07, 2.93e-07, 3.91e-07, 4.86e-07, 1.06e-06, 
                                                    2.93e-06, 3.48e-06, 2.89e-06, 9.89e-06, 1.38e-04, 
                                                    8.55e-05, 1.15e-04, 1.10e-03, 2.89e-04, 9.56e-04, 
                                                    3.17e-03, 6.70e-03, 1.90e-02, 5.90e-02, 1.15e-01, 
                                                    1.85e-01, 2.68e-01, 2.98e-01, 2.72e-01, 2.40e-01, 
                                                    1.92e-01, 1.35e-01, 9.24e-02, 6.10e-02, 3.68e-02, 
                                                    2.61e-02, 1.95e-02, 1.32e-02, 9.40e-03, 5.15e-03, 
                                                    3.60e-03, 3.40e-03, 3.80e-03, 4.60e-03 };

        // Water scattering and attenuation coefficient data (Morel 1988)
        std::vector<ScalarFloat> attn_wavelengths = {   400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 
                                                        450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 
                                                        500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 
                                                        550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 
                                                        600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 
                                                        650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 
                                                        700 };
        std::vector<ScalarFloat> attn_k = {     0.0209, 0.0200, 0.0196, 0.0189, 0.0183, 0.0182, 0.0171, 0.0170, 0.0168, 0.0166,
                                                0.0168, 0.0170, 0.0173, 0.0174, 0.0175, 0.0184, 0.0194, 0.0203, 0.0217, 0.0240,
                                                0.0271, 0.0320, 0.0384, 0.0445, 0.0490, 0.0505, 0.0518, 0.0543, 0.0568, 0.0615,
                                                0.0640, 0.0640, 0.0717, 0.0762, 0.0807, 0.0940, 0.1070, 0.1280, 0.1570, 0.2000,
                                                0.2530, 0.2790, 0.2960, 0.3030, 0.3100, 0.3150, 0.3200, 0.3250, 0.3300, 0.3400,
                                                0.3500, 0.3700, 0.4050, 0.4180, 0.4300, 0.4400, 0.4500, 0.4700, 0.5000, 0.5500,
                                                0.6500 };
        std::vector<ScalarFloat> attn_chi = {   0.1100, 0.1110, 0.1125, 0.1135, 0.1126, 0.1104, 0.1078, 0.1065, 0.1041, 0.0996,
                                                0.0971, 0.0939, 0.0896, 0.0859, 0.0823, 0.0788, 0.0746, 0.0726, 0.0690, 0.0660,
                                                0.0636, 0.0600, 0.0578, 0.0540, 0.0498, 0.0475, 0.0467, 0.0450, 0.0440, 0.0426,
                                                0.0410, 0.0400, 0.0390, 0.0375, 0.0360, 0.0340, 0.0330, 0.0328, 0.0325, 0.0330,
                                                0.0340, 0.0350, 0.0360, 0.0375, 0.0385, 0.0400, 0.0420, 0.0430, 0.0440, 0.0445,
                                                0.0450, 0.0460, 0.0475, 0.0490, 0.0515, 0.0520, 0.0505, 0.0440, 0.0390, 0.0340,
                                                0.0300 };
        std::vector<ScalarFloat> attn_e = {     0.668, 0.672, 0.680, 0.687, 0.693, 0.701, 0.707, 0.708, 0.707, 0.704,
                                                0.701, 0.699, 0.700, 0.703, 0.703, 0.703, 0.703, 0.704, 0.702, 0.700,
                                                0.700, 0.695, 0.690, 0.685, 0.680, 0.675, 0.670, 0.665, 0.660, 0.655,
                                                0.650, 0.645, 0.640, 0.630, 0.623, 0.615, 0.610, 0.614, 0.618, 0.622,
                                                0.626, 0.630, 0.634, 0.638, 0.642, 0.647, 0.653, 0.658, 0.663, 0.667,
                                                0.672, 0.677, 0.682, 0.687, 0.695, 0.697, 0.693, 0.665, 0.640, 0.620,
                                                0.600 };
        // IMPORTANT: This table uses the values provided by Morel, which are different than the ones from 6SV
        std::vector<ScalarFloat> molecular_scatter_coeffs = {   0.00618095, 0.00578095, 0.00547619, 0.00517619, 0.00492222, 
                                                                0.0046746 , 0.00447143, 0.00426825, 0.00406508, 0.0038619 , 
                                                                0.00365873, 0.00346667, 0.00331429, 0.0031619 , 0.00300952, 
                                                                0.00287143, 0.00276984, 0.00265238, 0.0025    , 0.00236508, 
                                                                0.00226349, 0.0021619 , 0.00206032, 0.00195873, 0.00185714, 
                                                                0.00177778, 0.00172698, 0.00167619, 0.0016254 , 0.0015746 , 
                                                                0.00152381, 0.00144603, 0.00134444, 0.0013    , 0.0013    , 
                                                                0.00126984, 0.00121905, 0.00116825, 0.00111746, 0.00107   , 
                                                                0.00102429, 0.00098556, 0.00095   , 0.0009181 , 0.00088762, 
                                                                0.00085714, 0.00082667, 0.00079619, 0.00076571, 0.00073937, 
                                                                0.00071397, 0.00069286, 0.00067254, 0.00065222, 0.0006319 , 
                                                                0.00061159, 0.00059127, 0.00057095, 0.00055063, 0.00053524, 
                                                                0.00052 };

        std::vector<ScalarFloat> molecular_scatter_coeffs_6s = {    0.0076, 0.0072, 0.0068, 0.0064, 0.0061,
                                                                    0.0058, 0.0055, 0.0052, 0.0049, 0.0047,
                                                                    0.0045, 0.0043, 0.0041, 0.0039, 0.0037,
                                                                    0.0036, 0.0034, 0.0033, 0.0031, 0.0030,
                                                                    0.0029, 0.0027, 0.0026, 0.0025, 0.0024,
                                                                    0.0023, 0.0022, 0.0022, 0.0021, 0.0020,
                                                                    0.0019, 0.0018, 0.0018, 0.0017, 0.0017,
                                                                    0.0016, 0.0016, 0.0015, 0.0015, 0.0014,
                                                                    0.0014, 0.0013, 0.0013, 0.0012, 0.0012,
                                                                    0.0011, 0.0011, 0.0010, 0.0010, 0.0010,
                                                                    0.0010, 0.0009, 0.0008, 0.0008, 0.0008,
                                                                    0.0007, 0.0007, 0.0007, 0.0007, 0.0007,
                                                                    0.0007 };

        std::vector<ScalarFloat> upwelling_transmittance_x = {  0.007578124952223935, 5.059528519788923, 10.013657992469575, 15.005149724840436, 20.024686390257536, 
                                                                25.00098028843596, 30.058706007033862, 35.03919747216452, 40.0728238728059, 45.04713640196841, 
                                                                50.039794940537526, 55.04214017207334, 60.063398916240075, 65.06474612486417, 70.04005667693839, 
                                                                75.01690340153361, 80.06547079223573, 85.03637830205237, 88.8382585234031 };
        std::vector<ScalarFloat> upwelling_transmittance_y = {  0.9308556978035891, 0.9277408752589185, 0.9222447960969848, 0.9147882496572222, 0.9062573374856766, 
                                                                0.8976495372171885, 0.8890513169464371, 0.881106879830576, 0.874009383435481, 0.8680223810801799, 
                                                                0.8632014091832066, 0.8602915434534224, 0.8586984105060734, 0.8588850205914815, 0.8605153432278858, 
                                                                0.8626100545315526, 0.8662813991317904, 0.8727440104682497, 0.8786882818983944 };

        std::vector<ScalarFloat> downwelling_transmittance_x = {    0.025306818275454522, 5.084330509589584, 10.08454946858249, 15.08177562124021, 19.040803613490898, 
                                                                    25.005251470430057, 30.029586482432702, 35.003043786737216, 40.089352384827755, 45.02512575838238, 
                                                                    50.021170926187246, 55.09609512370865, 60.09015486770292, 65.01502308091857, 70.09880344228476, 
                                                                    75.06536555552293, 80.0779836329192, 85.01990888616281, 88.93081235693548 };
        std::vector<ScalarFloat> downwelling_transmittance_y = {    0.9766129527157106, 0.9766016598604734, 0.9763386932360046, 0.9759268756311694, 0.9753692521943176, 
                                                                    0.9744350075622045, 0.9734044983929273, 0.9721974463179486, 0.97064054626727, 0.9687443754002456, 
                                                                    0.9665641214393463, 0.9639756558797935, 0.9609870543861947, 0.9573905793466994, 0.9531669207358379, 
                                                                    0.948072872705384, 0.9422679803838081, 0.9357601370709014, 0.9295328947161291 };

        // Initialize the vectors with azimuth/zenith points/weights
        m_azim_pts = dr::Array<ScalarFloat, 48>(  3.86099459e-03f, 2.03255633e-02f, 4.98740911e-02f, 9.23892368e-02f,
                                            1.47693486e-01f, 2.15555068e-01f, 2.95689413e-01f, 3.87760434e-01f,
                                            4.91381968e-01f, 6.06119396e-01f, 7.31491473e-01f, 8.66972347e-01f,
                                            1.01199377e+00f, 1.16594746e+00f, 1.32818769e+00f, 1.49803398e+00f,
                                            1.67477392e+00f, 1.85766620e+00f, 2.04594372e+00f, 2.23881677e+00f,
                                            2.43547638e+00f, 2.63509768e+00f, 2.83684340e+00f, 3.03986735e+00f,
                                            3.24331796e+00f, 3.44634190e+00f, 3.64808762e+00f, 3.84770893e+00f,
                                            4.04436853e+00f, 4.23724158e+00f, 4.42551910e+00f, 4.60841139e+00f,
                                            4.78515133e+00f, 4.95499761e+00f, 5.11723785e+00f, 5.27119154e+00f,
                                            5.41621296e+00f, 5.55169383e+00f, 5.67706591e+00f, 5.79180334e+00f,
                                            5.89542487e+00f, 5.98749589e+00f, 6.06763024e+00f, 6.13549182e+00f,
                                            6.19079607e+00f, 6.23331122e+00f, 6.26285974e+00f, 6.27932431e+00f);

        m_zenith_pts = dr::Array<ScalarFloat, 24>(0.00038299f, 0.00201104f, 0.00491196f, 0.00903877f, 
                                            0.01432379f, 0.02068026f, 0.02800382f, 0.03617421f, 
                                            0.04505728f, 0.05450717f, 0.06436872f, 0.07447999f, 
                                            0.08467496f, 0.09478623f, 0.10464777f, 0.11409766f, 
                                            0.12298073f, 0.13115113f, 0.13847468f, 0.14483116f,
                                            0.15011618f, 0.15424299f, 0.15714391f, 0.15877195f);

        m_azim_weights = dr::Array<ScalarFloat, 48>(  0.00315335f, 0.00732755f, 0.01147723f, 0.01557932f, 0.01961616f,
                                                0.02357076f, 0.02742651f, 0.03116723f, 0.03477722f, 0.03824135f,
                                                0.04154508f, 0.04467456f, 0.04761666f, 0.05035904f, 0.05289019f,
                                                0.0551995f , 0.05727729f, 0.05911484f, 0.06070444f, 0.06203942f,
                                                0.06311419f, 0.06392424f, 0.06446616f, 0.0647377f , 0.0647377f ,
                                                0.06446616f, 0.06392424f, 0.06311419f, 0.06203942f, 0.06070444f,
                                                0.05911484f, 0.05727729f, 0.0551995f , 0.05289019f, 0.05035904f,
                                                0.04761666f, 0.04467456f, 0.04154508f, 0.03824135f, 0.03477722f,
                                                0.03116723f, 0.02742651f, 0.02357076f, 0.01961616f, 0.01557932f,
                                                0.01147723f, 0.00732755f, 0.00315335f);

        m_zenith_weights = dr::Array<ScalarFloat, 24>(0.01234123f, 0.02853139f, 0.04427744f, 0.05929858f, 0.07334648f,
                                                0.08619016f, 0.09761865f, 0.10744427f, 0.11550567f, 0.12167047f,
                                                0.12583746f, 0.1279382f , 0.1279382f , 0.12583746f, 0.12167047f,
                                                0.11550567f, 0.10744427f, 0.09761865f, 0.08619016f, 0.07334648f,
                                                0.05929858f, 0.04427744f, 0.02853139f, 0.01234123f);
 
        // Construct distributions from the provided data sets
        m_effective_reflectance = IrregularContinuousDistribution<ScalarFloat>(
            wc_wavelengths.data(), wc_data.data(), wc_data.size()
        );

        m_ior_real = IrregularContinuousDistribution<ScalarFloat>(
            ior_wavelengths.data(), ior_real_data.data(), ior_real_data.size()
        );

        m_ior_imag = IrregularContinuousDistribution<ScalarFloat>(
            ior_wavelengths.data(), ior_cplx_data.data(), ior_cplx_data.size()
        );

        m_attn_k = IrregularContinuousDistribution<ScalarFloat>(
            attn_wavelengths.data(), attn_k.data(), attn_k.size()
        );

        m_attn_chi = IrregularContinuousDistribution<ScalarFloat>(
            attn_wavelengths.data(), attn_chi.data(), attn_chi.size()
        );

        m_attn_e = IrregularContinuousDistribution<ScalarFloat>(
            attn_wavelengths.data(), attn_e.data(), attn_e.size()
        );

        m_molecular_scatter_coeffs = IrregularContinuousDistribution<ScalarFloat>(
            attn_wavelengths.data(), molecular_scatter_coeffs.data(), molecular_scatter_coeffs.size()
        );

        m_molecular_scatter_coeffs_6s = IrregularContinuousDistribution<ScalarFloat>(
            attn_wavelengths.data(), molecular_scatter_coeffs_6s.data(), molecular_scatter_coeffs_6s.size()
        );

        m_upwelling_transmittance = IrregularContinuousDistribution<Float>(
            upwelling_transmittance_x.data(), upwelling_transmittance_y.data(), upwelling_transmittance_y.size()
        );

        m_downwelling_transmittance = IrregularContinuousDistribution<Float>(
            downwelling_transmittance_x.data(), downwelling_transmittance_y.data(), downwelling_transmittance_y.size()
        );

    }

    ScalarFloat effective_reflectance(const ScalarFloat &wavelength) const {
        return m_effective_reflectance.eval_pdf(wavelength);
    }

    ScalarFloat ior_real(const ScalarFloat &wavelength) const {
        return m_ior_real.eval_pdf(wavelength);
    }

    ScalarFloat ior_cplx(const ScalarFloat &wavelength) const {
        return m_ior_imag.eval_pdf(wavelength);
    }

    ScalarFloat attn_k(const ScalarFloat &wavelength) const {
        return m_attn_k.eval_pdf(wavelength);
    }

    ScalarFloat attn_chi(const ScalarFloat &wavelength) const {
        return m_attn_chi.eval_pdf(wavelength);
    }

    ScalarFloat attn_e(const ScalarFloat &wavelength) const {
        return m_attn_e.eval_pdf(wavelength);
    }

    ScalarFloat molecular_scatter_coeff(const ScalarFloat &wavelength) const {
        return m_molecular_scatter_coeffs.eval_pdf(wavelength);
    }

    ScalarFloat molecular_scatter_coeff_6s(const ScalarFloat &wavelength) const {
        return m_molecular_scatter_coeffs_6s.eval_pdf(wavelength);
    }

    Float upwelling_transmittance(const Float &wavelength) const {
        return m_upwelling_transmittance.eval_pdf(wavelength);
    }

    Float downwelling_transmittance(const Float &wavelength) const {
        return m_downwelling_transmittance.eval_pdf(wavelength);
    }

    Value azimuth_point(const Int32 &idx) const {
        return dr::gather<Value>(m_azim_pts, idx);
    }

    Value zenith_point(const Int32 &idx) const {
        return dr::gather<Value>(m_zenith_pts, idx);
    }

    Value azimuth_weight(const Int32 &idx) const {
        return dr::gather<Value>(m_azim_weights, idx);
    }

    Value zenith_weight(const Int32 &idx) const {
        return dr::gather<Value>(m_zenith_weights, idx);
    }

    AzimuthStorage azimuth_points() const {
        return m_azim_pts;
    }

    ZenithStorage zenith_points() const {
        return m_zenith_pts;
    }

    AzimuthStorage azimuth_weights() const {
        return m_azim_weights;
    }

    ZenithStorage zenith_weights() const {
        return m_zenith_weights;
    }

private:
    // Effective reflectance of whitecaps
    IrregularContinuousDistribution<ScalarFloat> m_effective_reflectance;

    // Real/Complex IOR of water (Hale & Querry 1973)
    IrregularContinuousDistribution<ScalarFloat> m_ior_real;
    IrregularContinuousDistribution<ScalarFloat> m_ior_imag;

    // Water scattering and attenuation coefficients (Morel 1988)
    IrregularContinuousDistribution<ScalarFloat> m_attn_k;
    IrregularContinuousDistribution<ScalarFloat> m_attn_chi;
    IrregularContinuousDistribution<ScalarFloat> m_attn_e;
    IrregularContinuousDistribution<ScalarFloat> m_molecular_scatter_coeffs;
    IrregularContinuousDistribution<ScalarFloat> m_molecular_scatter_coeffs_6s;

    // Upwelling and downwelling transmittance data
    IrregularContinuousDistribution<Float> m_upwelling_transmittance;
    IrregularContinuousDistribution<Float> m_downwelling_transmittance;

    // Azimuth/Zenith points and weights for quadrature
    AzimuthStorage m_azim_pts;
    ZenithStorage m_zenith_pts;
    
    AzimuthStorage m_azim_weights;
    ZenithStorage m_zenith_weights;
};

template<typename Float, typename Spectrum>
class OceanUtilities {
public:
    MI_IMPORT_TYPES()

    using AzimuthStorage = dr::Array<ScalarFloat, 48>;
    using ZenithStorage = dr::Array<ScalarFloat, 24>;
    using Value = dr::Array<ScalarFloat, 1>;

    OceanUtilities() : m_ocean_props() { }

    ScalarFloat eval_whitecap_coverage(const ScalarFloat &wind_speed) {
        return dr::clamp(m_monahan_alpha * dr::pow(wind_speed, m_monahan_lambda), 0.0f, 1.0f);
    }

    ScalarFloat eval_whitecaps(const ScalarFloat &wavelength, const ScalarFloat &wind_speed) {
        // Compute the fractional coverage of whitecaps
        ScalarFloat coverage = eval_whitecap_coverage(wind_speed);

        // Old
        ScalarFloat eff_reflectance = m_ocean_props.effective_reflectance(wavelength);

        // Compute the whitecap reflectance
        ScalarFloat whitecap_reflectance = coverage * eff_reflectance;

        return whitecap_reflectance;
    }

    Float eval_glint(const ScalarFloat &wavelength, 
                     const Vector3f &wi, const Vector3f &wo, 
                     const ScalarFloat &wind_direction, const ScalarFloat &wind_speed,
                     const ScalarFloat &chlorinity) {
        // Transform directions into azimuthal and zenithal angles
        Float theta_i = dr::acos(wi.z());
        Float theta_o = dr::acos(wo.z());
        Float phi_i = dr::atan2(wi.y(), wi.x());
        Float phi_o = dr::atan2(wo.y(), wo.x());

        return eval_sun_glint<Float>(wavelength, theta_i, theta_o, phi_i, phi_o, wind_direction, wind_speed, chlorinity);
    }

    Float eval_underlight(const ScalarFloat &wavelength, 
                          const Vector3f &wi, const Vector3f &wo,
                          const ScalarFloat &wind_direction, const ScalarFloat &wind_speed,
                          const ScalarFloat &chlorinity, const ScalarFloat &pigmentation) {
        // Transform directions into azimuthal and zenithal angles
        Float theta_i = dr::acos(wi.z());
        Float theta_o = dr::acos(wo.z());

        return eval_underlight_angular(wavelength, theta_i, theta_o, wind_direction, wind_speed, chlorinity, pigmentation);
    }

private:
    // Ocean properties
    OceanProperties<Float, Spectrum> m_ocean_props;

    // Simple constants
    ScalarFloat m_pi = dr::Pi<ScalarFloat>;
    ScalarFloat m_pi_half = m_pi / 2.0f;
    ScalarFloat m_pi_three_half = (3.0f * m_pi) / 2.0f;

    // Whitecap parameters
    ScalarFloat m_f_eff_base = 0.4f;
    ScalarFloat m_monahan_alpha = 2.95e-06f;
    ScalarFloat m_monahan_lambda = 3.52f;

    // Cox-Munk distribution parameters
    ScalarFloat m_c_40 = 0.40f;
    ScalarFloat m_c_22 = 0.12f;
    ScalarFloat m_c_04 = 0.23f;

    // Underlight parameters
    ScalarFloat m_underlight_alpha = 0.485f;

    template <typename Value>
    Value eval_cox_munk(const Value &phi_w, 
                        const Value &z_x, const Value &z_y,
                        const ScalarFloat &wind_speed) const {
        ScalarFloat sigma_c = 0.003f + 0.00192f * wind_speed;
        ScalarFloat sigma_u = 0.00316f * wind_speed;

        ScalarFloat c_21 = 0.01f - 0.0086f * wind_speed;
        ScalarFloat c_03 = 0.04f - 0.033f * wind_speed;

        Value xe = (dr::cos(phi_w) * z_x + dr::sin(phi_w) * z_y) / dr::sqrt(sigma_c);
        Value xn = (-dr::sin(phi_w) * z_x + dr::cos(phi_w) * z_y) / dr::sqrt(sigma_u);

        Value xe2 = xe * xe;
        Value xn2 = xn * xn;
        
        Value coef = 1.0f - (c_21 / 2.0f) * (xe2 - 1.0f) * xn - (c_03 / 6.0f) * (xn2 - 3.0f) * xn;
        coef = coef + (m_c_40 / 24.0f) * (xe2 * xe2 - 6.0f * xe2 + 3.0f);
        coef = coef + (m_c_04 / 24.0f) * (xn2 * xn2 - 6.0f * xn2 + 3.0f);
        coef = coef + (m_c_22 / 4.0f) * (xe2 - 1.0f) * (xn2 - 1.0f);
        
        Value prob = coef / 2.0f / dr::Pi<Value> / dr::sqrt(sigma_u) / dr::sqrt(sigma_c) * dr::exp(-(xe2 + xn2) / 2.0f);
        return prob;
    }

    template <typename Value>
    Value eval_fresnel(const ScalarFloat &n_real, const ScalarFloat &n_imag,
                       const Float &coschi, const Float &sinchi) const {
        Value s = (n_real * n_real) - (n_imag * n_imag) - (sinchi * sinchi);
        
        Value a_1 = dr::abs(s);
        Value a_2 = dr::sqrt(dr::sqr(s) + 4.0f * n_real * n_real * n_imag * n_imag);

        Value u = dr::sqrt(0.5f * dr::abs(a_1 + a_2));
        Value v = dr::sqrt(0.5f * dr::abs(a_2 - a_1));

        Value b_1 = (n_real * n_real - n_imag * n_imag) * coschi;
        Value b_2 = 2 * n_real * n_imag * coschi;

        Value right_squared = (dr::sqr(coschi - u) + v * v) / (dr::sqr(coschi + u) + v * v);
        Value left_squared = (dr::sqr(b_1 - u) + dr::sqr(b_2 + v)) / (dr::sqr(b_1 + u) + dr::sqr(b_2 - v));
        Value R = (right_squared + left_squared) * 0.5f;

        return R;
    }

    template <typename Value>
    Value eval_sun_glint(const ScalarFloat &wavelength, 
                         const Value &theta_i, const Value &theta_o,
                         const Value &phi_i, const Value &phi_o,
                         const ScalarFloat &wind_direction, const ScalarFloat &wind_speed,
                         const ScalarFloat &chlorinity) {
        // Implementation analog to 6SV
        Value phi = phi_i - phi_o;
        Value phi_w = phi_i - wind_direction;

        Value c_i = dr::cos(theta_i);
        Value c_o = dr::cos(theta_o);
        Value s_i = dr::sin(theta_i);
        Value s_o = dr::sin(theta_o);

        Value z_x = (-s_o * dr::sin(phi)) / (c_i + c_o);
        Value z_y = (s_i + s_o * dr::cos(phi)) / (c_i + c_o);

        // Tilt angle (rad)
        Value tan_tilt = dr::sqrt(z_x * z_x + z_y * z_y);
        Value tilt = dr::atan(tan_tilt);

        // Cox-Munk specular probability
        Value specular_prob = eval_cox_munk<Value>(phi_w, z_x, z_y, wind_speed);
        auto mask = Mask(specular_prob < 0.0f);
        specular_prob = dr::select(mask, 0.0f, specular_prob);

        Value cos_2_chi = c_o * c_i + s_o * s_i * dr::cos(phi);
        auto ge_1 = Mask(cos_2_chi > 1.0f);
        auto le_1 = Mask(cos_2_chi < -1.0f);
        
        cos_2_chi = dr::select(ge_1, 0.999999999f, cos_2_chi);
        cos_2_chi = dr::select(le_1, -0.999999999f, cos_2_chi);
        
        Value coschi = dr::sqrt(0.5f * (1.0f + cos_2_chi));
        Value sinchi = dr::sqrt(0.5f * (1.0f - cos_2_chi));

        auto ge_coschi = Mask(coschi > 1.0f);
        auto le_coschi = Mask(coschi < -1.0f);
        auto ge_sinchi = Mask(sinchi > 1.0f);
        auto le_sinchi = Mask(sinchi < -1.0f);

        coschi = dr::select(ge_coschi, 0.999999999f, coschi);
        coschi = dr::select(le_coschi, -0.999999999f, coschi);
        sinchi = dr::select(ge_sinchi, 0.999999999f, sinchi);
        sinchi = dr::select(le_sinchi, -0.999999999f, sinchi);

        // Fresnel coefficient
        ScalarFloat n_real = m_ocean_props.ior_real(wavelength) + friedman_sverdrup(chlorinity);
        ScalarFloat n_imag = m_ocean_props.ior_cplx(wavelength);
        Value fresnel_coeff = eval_fresnel<Value>(n_real, n_imag, coschi, sinchi);
        
        // Sun glint reflectance
        Value num = m_pi * fresnel_coeff * specular_prob;
        Value denom = 4.0f * c_i * c_o * dr::pow(dr::cos(tilt), 4.0f);

        return num / denom;
    }

    Float eval_underlight_angular(const ScalarFloat &wavelength, 
                          const Float &theta_i, const Float &theta_o,
                          const ScalarFloat &wind_direction, const ScalarFloat &wind_speed,
                          const ScalarFloat &chlorinity, const ScalarFloat &pigmentation) {
        // Analogue to 6SV, we return 0.0 if the wavelength is outside the range of [0.4, 0.7]
        auto outside_range = Mask(wavelength < 0.4f || wavelength > 0.7f);

        // Get IOR of water
        ScalarFloat n_real = m_ocean_props.ior_real(wavelength) + friedman_sverdrup(chlorinity);
        ScalarFloat n_imag = m_ocean_props.ior_cplx(wavelength);

        // Compute r_omega
        ScalarFloat r_om = r_omega(wavelength, pigmentation);

        // Upwelling and downwelling transmittance
        Float t_u = upwelling_transmittance_lut(theta_o);
        Float t_d = downwelling_transmittance_lut(theta_i);

        downwelling_transmittance_quadrature(wavelength, theta_i, wind_direction, wind_speed, chlorinity);

        // Compute the underlight term
        Float underlight = (1.0f / (dr::sqr(n_real) + dr::sqr(n_imag))) * (r_om * t_u * t_d) / (1.0f - m_underlight_alpha * r_om);

        return dr::select(outside_range, 0.0f, underlight);
    }

    // Correction to the IOR of water according to Friedman (1969) and Sverdrup (1942)
    ScalarFloat friedman_sverdrup(const ScalarFloat &chlorinity) {
        return 0.00017492711f * (0.03f + 1.805f * chlorinity);
    }

    Float downwelling_transmittance_quadrature(const ScalarFloat &wavelength, 
                                               const ScalarFloat &theta_i,
                                               const ScalarFloat &wind_direction, const ScalarFloat &wind_speed,
                                               const ScalarFloat &chlorinity) {  
        // Construct azimuth storage for wavelength and theta_i
        AzimuthStorage theta_is = dr::Array<ScalarFloat, 48>(theta_i);
        AzimuthStorage phi_is = dr::Array<ScalarFloat, 48>(0.0f);

        // Perform quadrature loop
        for (int zenith_idx = 0; zenith_idx < 24; zenith_idx++) {
            Value zenith = m_ocean_props.zenith_point(Int32(zenith_idx));
            AzimuthStorage azimuths = m_ocean_props.azimuth_points();

            Value zenith_weight = m_ocean_props.zenith_weight(Int32(zenith_idx));
            AzimuthStorage azimuth_weights = m_ocean_props.azimuth_weights();

            // Access the zenith point and weight 
            ScalarFloat zenith_f = zenith[0];
            ScalarFloat zenith_weight_f = zenith_weight[0];
            
            // To obtain the final weight, we multiply the zenith weight with all azimuthal weights
            AzimuthStorage final_weights = zenith_weight_f * azimuth_weights;

            // Finally, we need an array of zeniths of length 48
            AzimuthStorage zeniths = dr::Array<ScalarFloat, 48>(zenith_f);

            // Evaluate the glint term for all azimuthal points and the current zenith point
            AzimuthStorage glint = eval_sun_glint<AzimuthStorage>(wavelength, theta_is, zeniths, phi_is, azimuths, wind_direction, wind_speed, chlorinity);

            Log(Warn, "Glint: %s", glint);

            // Multiply by the weights
            //AzimuthStorage weighted_glint = final_weights * glint;

            //Log(Warn, "Weighted glint: %s", weighted_glint);
        }
    
        return 0.0f;
    }

    Float upwelling_transmittance_quadrature(const Float &theta_o) {
        
    }

    Float downwelling_transmittance_polynomial(const Float &theta_i) {
        // Evaluate fitted downwelling polynomial
        // Coefficients: 0.00882913 -0.03803999  0.02232457 -0.00533469  0.97575765
        return 0.00882913 * dr::pow(theta_i, 4.0f)
                - 0.03803999 * dr::pow(theta_i, 3.0f)
                + 0.02232457 * dr::pow(theta_i, 2.0f)
                - 0.00533469 * theta_i
                + 0.97575765;
    }

    Float upwelling_transmittance_polynomial(const Float &theta_o) {
        // Evaluate fitted downwelling polynomial
        // Coefficients: -0.03694741  0.17051226 -0.18538372 -0.01917329  0.93060225
        return -0.03694741 * dr::pow(theta_o, 4.0f) 
            + 0.17051226 * dr::pow(theta_o, 3.0f)
            - 0.18538372 * dr::pow(theta_o, 2.0f)
            - 0.01917329 * theta_o
            + 0.93060225; 
    }

    Float downwelling_transmittance_lut(const Float &theta_i) {
        return m_ocean_props.downwelling_transmittance(theta_i);
    }

    Float upwelling_transmittance_lut(const Float &theta_o) {
        return m_ocean_props.upwelling_transmittance(theta_o);
    }

    ScalarFloat r_omega(const ScalarFloat &wavelength,
                        const ScalarFloat &pigmentation) {
        ScalarFloat wavelength_nm = wavelength * 1000.0f;
        ScalarFloat pigment_log = dr::log(pigmentation) / dr::log(10.0f);

        // Backscattering coefficient
        ScalarFloat molecular_scatter_coeff = m_ocean_props.molecular_scatter_coeff_6s(wavelength_nm);
        ScalarFloat scattering_coeff = 0.30f * dr::pow(pigmentation, 0.62);
        ScalarFloat backscatter_ratio = 0.002f + 0.02f * (0.5f - 0.25f * pigment_log) * (550.0 / wavelength_nm);
        ScalarFloat backscatter_coeff = 0.5f * molecular_scatter_coeff + scattering_coeff * backscatter_ratio;

        // (Diffuse) attenuation coefficient
        ScalarFloat k = m_ocean_props.attn_k(wavelength_nm);
        ScalarFloat chi = m_ocean_props.attn_chi(wavelength_nm);
        ScalarFloat e = m_ocean_props.attn_e(wavelength_nm);
        ScalarFloat attn_coeff = k + chi * dr::pow(pigmentation, e);

        // If any of the coefficients is zero, we return zero
        if (backscatter_coeff == 0.0f || attn_coeff == 0.0f)
            return 0.0f;

        // Iterative computation of the reflectance
        ScalarFloat u = 0.75f;
        ScalarFloat r_omega = 0.33f * backscatter_coeff / u / attn_coeff;

        bool converged = false;
        while (!converged) {
            // Update u
            u = (0.9f * (1.0f - r_omega)) / (1.0f + 2.25f * r_omega);

            // Update reflectance
            ScalarFloat r_omega_new = 0.33f * backscatter_coeff / (u * attn_coeff);

            // Create a mask that marks the converged values
            if (dr::abs((r_omega_new - r_omega) / r_omega_new) < 0.0001f) {
                converged = true;
                break;
            }

            // Update reflectance ONLY for non-converged values
            r_omega = r_omega_new;
        }

        return r_omega;
    }

};

#endif // OCEAN_PROPS

template <typename Float, typename Spectrum>
class OceanicBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    OceanicBSDF(const Properties &props) : Base(props) {
        // Retrieve the parameters used in 6SV
        m_channel = props.get<ScalarInt32>("channel");
        m_type = props.get<ScalarInt32>("visual_type");
        m_wavelength = props.get<ScalarFloat>("wavelength");
        m_wind_speed = props.get<ScalarFloat>("wind_speed");
        m_wind_direction = props.get<ScalarFloat>("wind_direction");
        m_chlorinity = props.get<ScalarFloat>("chlorinity");
        m_pigmentation = props.get<ScalarFloat>("pigmentation");
        m_shininess = props.get<ScalarFloat>("shininess");

        // Initialize the ocean utilities
        m_ocean_utils = new OceanUtilities<Float, Spectrum>();

        // Set the BSDF flags
        // => Whitecap and underlight reflectance is "diffuse"
        m_components.push_back(BSDFFlags::DiffuseReflection | 
                               BSDFFlags::FrontSide);
    
        // => Sun glint reflectance at the water surface is "specular"
        m_components.push_back(BSDFFlags::GlossyReflection | 
                               BSDFFlags::FrontSide | BSDFFlags::BackSide);

        // Set all the flags
        for (auto c : m_components)
            m_flags |= c;
        dr::set_attr(this, "flags", m_flags);
    }

    Float eval_whitecaps() const {
        return m_ocean_utils->eval_whitecaps(m_wavelength, m_wind_speed);
    }

    Float eval_glint(const Vector3f &wi, const Vector3f &wo) const {
        return m_ocean_utils->eval_glint(m_wavelength, wi, wo, m_wind_direction, m_wind_speed, m_chlorinity);
    }

    Float eval_underlight(const Vector3f &wi, const Vector3f &wo) const {
        return m_ocean_utils->eval_underlight(m_wavelength, wi, wo, m_wind_direction, m_wind_speed, m_chlorinity, m_pigmentation);
    }

    Float eval_ocean(const Vector3f &wi, const Vector3f &wo) const {
        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float whitecap_reflectance = eval_whitecaps();
        Float glint_reflectance = eval_glint(wi, wo);
        Float underlight_reflectance = eval_underlight(wi, wo);

        return (coverage * whitecap_reflectance) 
            + (1 - coverage) * glint_reflectance
            + (1 - (coverage * whitecap_reflectance)) * underlight_reflectance;
    }

    Float eval_blinn_phong(const Vector3f &wi, const Vector3f &wo, const Normal3f &normal) const {
        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float factor = (m_shininess + 2.0f) / (2.0f * dr::Pi<Float>);
        Vector3f half = dr::normalize(wi + wo);
        Float dot = dr::dot(half, normal);

        // Blinn-phong => clamp dot to above zero
        dot = dr::clamp(dot, 0.0f, 1.0f);

        Float phong = dr::pow(dot, m_shininess);

        return coverage + (1 - coverage) * factor * phong;
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active)) || (!has_diffuse && !has_specular))
            return { bs, 0.f };

        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float prob_diff = coverage;

        Mask sample_diffuse = active & (sample1 < prob_diff),
             sample_glint = active && !sample_diffuse;

        if (dr::any_or<true>(sample_diffuse)) {
            // In the case of sampling the diffuse component, the outgoing direction
            // is sampled from a cosine-weighted hemisphere.
            dr::masked(bs.wo, sample_diffuse) = warp::square_to_cosine_hemisphere(sample2);
            dr::masked(bs.sampled_component, sample_diffuse) = 0;
            dr::masked(bs.sampled_type, sample_diffuse) = +BSDFFlags::DiffuseReflection;
        }

        if (dr::any_or<true>(sample_glint)) {
            // For Blinn-Phong, we need to sample the half-vector
            Float ksi_1 = sample2.x(),
                  ksi_2 = sample2.y();

            Float phi_h = dr::TwoPi<Float> * ksi_1;
            Float theta_h = dr::acos(dr::pow(ksi_2, 1.f / (m_shininess + 2.f)));

            Vector3f half = dr::normalize(Vector3f(dr::sin(theta_h) * dr::cos(phi_h),
                                                 dr::sin(theta_h) * dr::sin(phi_h),
                                                 dr::cos(theta_h)));

            Vector3f wo = 2.f * dr::dot(si.wi, half) * half - si.wi;

            // In the case of sampling the glint component, the outgoing direction
            // is sampled using the Blinn-Phong distribution.
            dr::masked(bs.wo, sample_glint) = wo; 
            dr::masked(bs.sampled_component, sample_glint) = 1;
            dr::masked(bs.sampled_type, sample_glint) = +BSDFFlags::GlossyReflection;
        }

        bs.pdf = pdf(ctx, si, bs.wo, active);
        bs.eta = 1.f;

        UnpolarizedSpectrum value = eval_ocean(si.wi, bs.wo) * Frame3f::cos_theta(bs.wo) / bs.pdf;

        return { bs, (depolarizer<Spectrum>(value)) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        bool has_underlight = ctx.is_enabled(BSDFFlags::DiffuseTransmission, 2);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap && !has_glint && !has_underlight))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);
        UnpolarizedSpectrum blinn(0.f);
        UnpolarizedSpectrum whitecap_reflectance(0.f);
        UnpolarizedSpectrum glint_reflectance(0.f);
        UnpolarizedSpectrum underlight_reflectance(0.f);
        
        // Get the reflected directions
        auto is_reflect = Mask(dr::eq(dr::sign(cos_theta_i), dr::sign(cos_theta_o))) && active;

        if (has_whitecap)
            // If whitecaps are enabled, compute the whitecap reflectance 
            whitecap_reflectance = eval_whitecaps();
        
        if (has_glint)
            // If sun glint is enabled, compute the glint reflectance
            glint_reflectance = eval_glint(si.wi, wo);
    
        if (has_underlight)
            underlight_reflectance = eval_underlight(si.wi, wo);
        
        // Combine the results
        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);

        // For debugging purposes, the channel indicates what term of the BRDF to evaluate
        switch (m_channel)
        {
            case 0:
                result[is_reflect] = (whitecap_reflectance) & active;
                break;
            case 1:
                result[is_reflect] = ((1 - coverage) * glint_reflectance) & active;
                break;
            case 2:
                result[is_reflect] = ((1 - (coverage * whitecap_reflectance)) * underlight_reflectance) & active;
                break;        
            default:
                result[is_reflect] = whitecap_reflectance
                    + (1 - coverage) * glint_reflectance
                    + (1 - (coverage * whitecap_reflectance)) * underlight_reflectance;
                break;
        }

        // Cosine foreshortening factor
        result[is_reflect] *= cos_theta_o;

        // Compute the Blinn-Phong term
        blinn[is_reflect] = eval_blinn_phong(si.wi, wo, si.n);

        // Compute normalized values over the vector
        auto normalized_r = result / dr::max(dr::max(result));
        auto normalized_b = blinn / dr::max(dr::max(blinn));

        // For debugging and visualization purposes, we support multiple visualisation types
        switch (m_type) 
        {
            // Visualize the ocean reflectance
            case 0:
                break;

            // Visualize the Blinn-Phong term
            case 1:
                result = blinn;
                break;

            // Visualize the normalized ocean reflectance
            case 2:
                result = normalized_r;
                break;

            // Visualize the difference between the normalized ocean and 
            // Blinn-Phong reflectance
            case 3:
                result = normalized_r - normalized_b;
                break;
        
            default:
                break;
        }

        return depolarizer<Spectrum>(result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        if (unlikely((!has_diffuse && !has_specular) || dr::none_or<false>(active)))
            return 0.f;

        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float weight_diffuse = coverage,
              weight_specular = (1 - coverage);

        // Check if the normal has only zeros. If this is the case, use a default normal
        Vector3f normal = si.n;
        if (dr::all(normal == 0.f))
            normal = Vector3f(0.f, 0.f, 1.f);

        Vector3f half = dr::normalize(si.wi + wo);
        Float projection = dr::dot(half, normal);
        Float D = ((m_shininess + 2.0f) / dr::TwoPi<Float>) * dr::pow(projection, m_shininess);

        // We multiply the probability of the specular lobe with the pdf of 
        // the Blinn-Phong distribution and the probability of the diffuse lobe
        // with the pdf of the cosine-weighted hemisphere.
        Float pdf_diffuse = warp::square_to_cosine_hemisphere_pdf(wo),
              pdf_specular = (D * projection) / (4.0f * dr::dot(si.wi, half));

        Float pdf = weight_diffuse * pdf_diffuse + weight_specular * pdf_specular;

        // If the outgoing direction is in the lower hemisphere, we return zero
        Float cos_theta_o = Frame3f::cos_theta(wo);

        return dr::select(cos_theta_o > 0.f, pdf, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanicLegacy[" << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << std::endl
            << "  wind_direction = " << string::indent(m_wind_direction) << std::endl
            << "  chlorinity = " << string::indent(m_chlorinity) << std::endl
            << "  pigmentation = " << string::indent(m_pigmentation) << std::endl
            << "  shininess = " << string::indent(m_shininess) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    // User-provided fields
    ScalarInt32 m_channel;
    ScalarInt32 m_type;
    ScalarFloat m_wavelength;
    ScalarFloat m_wind_speed;
    ScalarFloat m_wind_direction;
    ScalarFloat m_chlorinity;
    ScalarFloat m_pigmentation;
    ScalarFloat m_shininess;

    // Fields used to compute whitecap reflectance
    OceanUtilities<Float, Spectrum> *m_ocean_utils;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)