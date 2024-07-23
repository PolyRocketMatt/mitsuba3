#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <drjit/dynamic.h>

NAMESPACE_BEGIN(mitsuba)

//using FloatX = dr::DynamicArray<float>;

// Header content - Potentially move somewhere else later
#ifndef OCEAN_PROPS
#define OCEAN_PROPS

template<typename Float, typename Spectrum>
class OceanProperties {
public:
    MI_IMPORT_TYPES()

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

        // Construct distributions from the provided data sets
        m_effective_reflectance = IrregularContinuousDistribution<Float>(
            wc_wavelengths.data(), wc_data.data(), wc_data.size()
        );

        m_ior_real = IrregularContinuousDistribution<Float>(
            ior_wavelengths.data(), ior_real_data.data(), ior_real_data.size()
        );

        m_ior_imag = IrregularContinuousDistribution<Float>(
            ior_wavelengths.data(), ior_cplx_data.data(), ior_cplx_data.size()
        );

        m_attn_k = IrregularContinuousDistribution<Float>(
            attn_wavelengths.data(), attn_k.data(), attn_k.size()
        );

        m_attn_chi = IrregularContinuousDistribution<Float>(
            attn_wavelengths.data(), attn_chi.data(), attn_chi.size()
        );

        m_attn_e = IrregularContinuousDistribution<Float>(
            attn_wavelengths.data(), attn_e.data(), attn_e.size()
        );

        m_molecular_scatter_coeffs = IrregularContinuousDistribution<Float>(
            attn_wavelengths.data(), molecular_scatter_coeffs.data(), molecular_scatter_coeffs.size()
        );
    }

    Float ior_real(const Float &wavelength) const {
        return m_ior_real.eval_pdf(wavelength);
    }

    Float ior_cplx(const Float &wavelength) const {
        return m_ior_imag.eval_pdf(wavelength);
    }

    Float effective_reflectance(const Float &wavelength) const {
        return m_effective_reflectance.eval_pdf(wavelength);
    }

    Float attn_k(const Float &wavelength) const {
        return m_attn_k.eval_pdf(wavelength);
    }

    Float attn_chi(const Float &wavelength) const {
        return m_attn_chi.eval_pdf(wavelength);
    }

    Float attn_e(const Float &wavelength) const {
        return m_attn_e.eval_pdf(wavelength);
    }

    Float molecular_scatter_coeff(const Float &wavelength) const {
        return m_molecular_scatter_coeffs.eval_pdf(wavelength);
    }

private:
    // Effective reflectance of whitecaps
    IrregularContinuousDistribution<Float> m_effective_reflectance;

    // Real/Complex IOR of water (Hale & Querry 1973)
    IrregularContinuousDistribution<Float> m_ior_real;
    IrregularContinuousDistribution<Float> m_ior_imag;

    // Water scattering and attenuation coefficients (Morel 1988)
    IrregularContinuousDistribution<Float> m_attn_k;
    IrregularContinuousDistribution<Float> m_attn_chi;
    IrregularContinuousDistribution<Float> m_attn_e;
    IrregularContinuousDistribution<Float> m_molecular_scatter_coeffs;
};

template<typename Float, typename Spectrum>
class OceanUtilities {
public:
    MI_IMPORT_TYPES()

    using FloatStorage = dr::DynamicArray<Float>;

    OceanUtilities() : m_ocean_props() { }

    Float eval_whitecap_coverage(const Float &wind_speed) {
        return dr::clamp(m_monahan_alpha * dr::pow(wind_speed, m_monahan_lambda), 0.0f, 1.0f);
    }

    Float eval_whitecaps(const Float &wavelength, const Float &wind_speed) {
        // Compute the fractional coverage of whitecaps
        Float coverage = eval_whitecap_coverage(wind_speed);

        // Compute the efficiency factor
        Float efficiency = m_f_eff_base;

        // Compute the effective reflectance
        Float eff_reflectance = m_ocean_props.effective_reflectance(wavelength) + 0.10f;

        // Compute the whitecap reflectance
        Float whitecap_reflectance = coverage * efficiency * eff_reflectance;

        return whitecap_reflectance;
    }

    Float eval_cox_munk(const Float &phi_w, 
                        const Float &z_x, const Float &z_y,
                        const Float &wind_speed) const {
        Float sigma_c = 0.003f + 0.00192f * wind_speed;
        Float sigma_u = 0.00316f * wind_speed;

        Float c_21 = 0.01f - 0.0086f * wind_speed;
        Float c_03 = 0.04f - 0.033f * wind_speed;

        Float xe = (safe_cos(phi_w) * z_x + safe_sin(phi_w) * z_y) / dr::sqrt(sigma_c);
        Float xn = (-safe_sin(phi_w) * z_x + safe_cos(phi_w) * z_y) / dr::sqrt(sigma_u);

        Float xe2 = xe * xe;
        Float xn2 = xn * xn;
        
        Float coef = 1.0f - (c_21 / 2.0f) * (xe2 - 1.0f) * xn - (c_03 / 6.0f) * (xn2 - 3.0f) * xn;
        coef = coef + (m_c_40 / 24.0f) * (xe2 * xe2 - 6.0f * xe2 + 3.0f);
        coef = coef + (m_c_04 / 24.0f) * (xn2 * xn2 - 6.0f * xn2 + 3.0f);
        coef = coef + (m_c_22 / 4.0f) * (xe2 - 1.0f) * (xn2 - 1.0f);
        
        Float prob = coef / 2.0f / dr::Pi<Float> / dr::sqrt(sigma_u) / dr::sqrt(sigma_c) * dr::exp(-(xe2 + xn2) / 2.0f);
        return prob;
    }

    Float eval_fresnel(const Float &n_real, const Float &n_imag,
                       const Float &coschi, const Float &sinchi) const {
        Float s = n_real * n_real - n_imag * n_imag - sinchi * sinchi;
        
        Float a_1 = dr::abs(s);
        Float a_2 = dr::sqrt(dr::sqr(s) + 4.0f * n_real * n_real * n_imag * n_imag);

        Float u = dr::sqrt(0.5f * dr::abs(a_1 + a_2));
        Float v = dr::sqrt(0.5f * dr::abs(a_2 - a_1));

        Float b_1 = (n_real * n_real - n_imag * n_imag) * coschi;
        Float b_2 = 2 * n_real * n_imag * coschi;

        Float right_squared = (dr::sqr(coschi - u) + v * v) / (dr::sqr(coschi + u) + v * v);
        Float left_squared = (dr::sqr(b_1 - u) + dr::sqr(b_2 + v)) / (dr::sqr(b_1 + u) + dr::sqr(b_2 - v));
        Float R = (right_squared + left_squared) / 2.0f;

        return R;
    }

    Float eval_glint(const Float &wavelength, 
                     const Vector3f &wi, const Vector3f &wo, 
                     const Float &wind_direction, const Float &wind_speed,
                     const Float &chlorinity) {
        // Transform directions into azimuthal and zenithal angles
        Float theta_i = dr::acos(wi.z());
        Float theta_o = dr::acos(wo.z());
        Float phi_i = dr::atan2(wi.y(), wi.x());
        Float phi_o = dr::atan2(wo.y(), wo.x());

        return eval_glint_internal(wavelength, theta_i, theta_o, phi_i, phi_o, wind_direction, wind_speed, chlorinity);
    }

    Float eval_glint_internal(const Float &wavelength, 
                     const Float &theta_i, const Float &theta_o,
                     const Float &phi_i, const Float &phi_o,
                     const Float &wind_direction, const Float &wind_speed,
                     const Float &chlorinity) {
        // Implementation analog to 6SV
        Float phi = phi_i - phi_o;
        Float phi_w = phi_i - wind_direction;

        // TODO: Make sure to only consider angles between ]0, pi/2[
        Float c_i = safe_cos(theta_i);
        Float c_o = safe_cos(m_pi_half - theta_o);
        Float s_i = safe_sin(theta_i);
        Float s_o = safe_sin(m_pi_half - theta_o);

        Float z_x = -s_o * safe_sin(phi) / (c_i + c_o);
        Float z_y = (s_i + s_o * safe_cos(phi)) / (c_i + c_o);

        // Tilt angle (rad)
        Float tan_tilt = dr::sqrt(z_x * z_x + z_y * z_y);
        Float tilt = dr::atan(tan_tilt);

        // Cox-Munk specular probability
        Float specular_prob = eval_cox_munk(phi_w, z_x, z_y, wind_speed);
        auto mask = Mask(specular_prob < 0.0f);
        specular_prob = dr::select(mask, 0.0f, specular_prob);

        Float cos_2_chi = c_o * c_i + s_o * s_i * dr::cos(phi);
        auto ge_1 = Mask(cos_2_chi > 1.0f);
        auto le_1 = Mask(cos_2_chi < -1.0f);
        
        cos_2_chi = dr::select(ge_1, 0.999999999f, cos_2_chi);
        cos_2_chi = dr::select(le_1, -0.999999999f, cos_2_chi);
        
        Float coschi = dr::sqrt(0.5f * (1.0f + cos_2_chi));
        Float sinchi = dr::sqrt(0.5f * (1.0f - cos_2_chi));

        auto ge_coschi = Mask(coschi > 1.0f);
        auto le_coschi = Mask(coschi < -1.0f);
        auto ge_sinchi = Mask(sinchi > 1.0f);
        auto le_sinchi = Mask(sinchi < -1.0f);

        coschi = dr::select(ge_coschi, 0.999999999f, coschi);
        coschi = dr::select(le_coschi, -0.999999999f, coschi);
        sinchi = dr::select(ge_sinchi, 0.999999999f, sinchi);
        sinchi = dr::select(le_sinchi, -0.999999999f, sinchi);

        // Fresnel coefficient
        Float n_real = m_ocean_props.ior_real(wavelength) + friedman_sverdrup(chlorinity);
        Float n_imag = m_ocean_props.ior_cplx(wavelength);
        Float fresnel_coeff = eval_fresnel(n_real, n_imag, coschi, sinchi);
        
        // Compute reflectance
        Float num = m_pi * specular_prob * fresnel_coeff;
        Float denom = 4.0f * c_i * c_o * dr::pow(dr::cos(tilt), 4.0f);

        return num / denom;
    }

    Float eval_underlight(const Float &wavelength, 
                          const Vector3f &wi, const Vector3f &wo,
                          const Float &wind_direction, const Float &wind_speed, 
                          const Float &chlorinity, const Float &pigmentation) {
        // Transform directions into azimuthal and zenithal angles
        Float theta_i = dr::acos(wi.z());
        Float theta_o = dr::acos(wo.z());
        Float phi_i = dr::atan2(wi.y(), wi.x());
        Float phi_o = dr::atan2(wo.y(), wo.x());

        return eval_underlight(wavelength, theta_i, theta_o, phi_i, phi_o, wind_direction, wind_speed, chlorinity, pigmentation);
    }

    Float eval_underlight(const Float &wavelength, 
                          const Float &theta_i, const Float &theta_o,
                          const Float &phi_i, const Float &phi_o,
                          const Float &wind_direction, const Float &wind_speed, 
                          const Float &chlorinity, const Float &pigmentation) {
        // Analogue to 6SV, we return 0.0 if the wavelength is outside the range of [0.4, 0.7]
        auto mask = Mask(wavelength >= 0.4f || wavelength <= 0.7f);

        // Get IOR of water
        Float n_real = m_ocean_props.ior_real(wavelength) + friedman_sverdrup(chlorinity);
        Float n_imag = m_ocean_props.ior_cplx(wavelength);

        // Compute r_omega
        Float r_om = r_omega(wavelength, pigmentation);

        // Upwelling and downwelling transmittance
        Float t_u = upwelling_transmittance_fast(theta_o);
        Float t_d = downwelling_transmittance_fast(theta_i);

        Log(Warn, "R_Omega: %f", r_om);
        Log(Warn, "T_U: %f", t_u);
        Log(Warn, "T_D: %f", t_d);

        // Compute the underlight term
        Float underlight = (1.0f / (dr::sqr(n_real) + dr::sqr(n_imag))) * (r_om * t_u * t_d) / (1.0f - m_underlight_alpha * r_om);

        return dr::select(mask, underlight, 0.0f);
    }

private:
    // Ocean properties
    OceanProperties<Float, Spectrum> m_ocean_props;

    // Simple constants
    ScalarFloat m_pi = dr::Pi<ScalarFloat>;
    ScalarFloat m_pi_half = m_pi / 2.0f;
    ScalarFloat m_pi_three_half = (3.0f * m_pi) / 2.0f;
    ScalarFloat m_trig_eps_cos = 0.01745154894888401f;
    ScalarFloat m_trig_eps_sin = 0.01745154894888401f;

    // Whitecap parameters
    ScalarFloat m_f_eff_base = 0.4f;
    ScalarFloat m_monahan_alpha = 2.951f * 1e-6;
    ScalarFloat m_monahan_lambda = 3.52f;

    // Cox-Munk distribution parameters
    ScalarFloat m_c_40 = 0.40f;
    ScalarFloat m_c_22 = 0.12f;
    ScalarFloat m_c_04 = 0.23f;

    // Underlight parameters
    ScalarFloat m_underlight_alpha = 0.485f;

    ScalarFloat azimuth_pts[48] = { 
       -0.99877101f, -0.99353017f, -0.98412458f, -0.97059159f, -0.9529877f ,
       -0.93138669f, -0.90587914f, -0.87657202f, -0.84358826f, -0.8070662f ,
       -0.76715903f, -0.72403413f, -0.67787238f, -0.6288674f , -0.57722473f,
       -0.52316097f, -0.4669029f , -0.40868648f, -0.34875589f, -0.28736249f,
       -0.22476379f, -0.16122236f, -0.0970047f , -0.03238017f,  0.03238017f,
        0.0970047f ,  0.16122236f,  0.22476379f,  0.28736249f,  0.34875589f,
        0.40868648f,  0.4669029f ,  0.52316097f,  0.57722473f,  0.6288674f ,
        0.67787238f,  0.72403413f,  0.76715903f,  0.8070662f ,  0.84358826f,
        0.87657202f,  0.90587914f,  0.93138669f,  0.9529877f ,  0.97059159f,
        0.98412458f,  0.99353017f,  0.99877101f };
    ScalarFloat azimuth_weights[48] = {
        0.00315335f, 0.00732755f, 0.01147723f, 0.01557932f, 0.01961616f,
        0.02357076f, 0.02742651f, 0.03116723f, 0.03477722f, 0.03824135f,
        0.04154508f, 0.04467456f, 0.04761666f, 0.05035904f, 0.05289019f,
        0.0551995f , 0.05727729f, 0.05911484f, 0.06070444f, 0.06203942f,
        0.06311419f, 0.06392424f, 0.06446616f, 0.0647377f , 0.0647377f ,
        0.06446616f, 0.06392424f, 0.06311419f, 0.06203942f, 0.06070444f,
        0.05911484f, 0.05727729f, 0.0551995f , 0.05289019f, 0.05035904f,
        0.04761666f, 0.04467456f, 0.04154508f, 0.03824135f, 0.03477722f,
        0.03116723f, 0.02742651f, 0.02357076f, 0.01961616f, 0.01557932f,
        0.01147723f, 0.00732755f, 0.00315335f
    };
    ScalarFloat zenith_pts[48] = {
       -0.99518722f, -0.97472856f, -0.93827455f, -0.88641553f, -0.82000199f,
       -0.74012419f, -0.64809365f, -0.54542147f, -0.43379351f, -0.31504268f,
       -0.19111887f, -0.06405689f,  0.06405689f,  0.19111887f,  0.31504268f,
        0.43379351f,  0.54542147f,  0.64809365f,  0.74012419f,  0.82000199f,
        0.88641553f,  0.93827455f,  0.97472856f,  0.99518722f
    };
    ScalarFloat zenith_weights[24] = {
        0.01234123f, 0.02853139f, 0.04427744f, 0.05929858f, 0.07334648f,
        0.08619016f, 0.09761865f, 0.10744427f, 0.11550567f, 0.12167047f,
        0.12583746f, 0.1279382f , 0.1279382f , 0.12583746f, 0.12167047f,
        0.11550567f, 0.10744427f, 0.09761865f, 0.08619016f, 0.07334648f,
        0.05929858f, 0.04427744f, 0.02853139f, 0.01234123f
    };
    

    // "Safe" version of the cosine function
    Float safe_cos(const Float &angle) const {
        // Fix less-than-zero angles
        auto test = Mask(angle < 0.0f);
        Float corrected_angle = dr::select(test, angle + 2.0f * m_pi, angle);

        // TODO: Angles > 2pi
        Float angle_cos = dr::cos(corrected_angle);

        // Vectorized Clamp
        auto mask = Mask(corrected_angle > m_pi_half && corrected_angle < m_pi_three_half);
        return dr::select(mask, 
            dr::clamp(angle_cos, -1.0f, -m_trig_eps_cos), 
            dr::clamp(angle_cos, m_trig_eps_cos, 1.0f));
    }

    // "Safe" version of the sine function
    Float safe_sin(const Float &angle) const {
        // Fix less-than-zero angles
        auto test = Mask(angle < 0.0f);
        Float corrected_angle = dr::select(test, angle + 2.0f * m_pi, angle);

        // TODO: Angles > 2pi
        Float angle_sine = dr::sin(corrected_angle);

        // Vectorized Clamp
        auto mask = Mask(corrected_angle > m_pi);
        return dr::select(mask, 
            dr::clamp(angle_sine, -1.0f, -m_trig_eps_sin),
            dr::clamp(angle_sine, m_trig_eps_sin, 1.0f));
    }

    // Correction to the IOR of water according to Friedman (1969) and Sverdrup (1942)
    Float friedman_sverdrup(const Float &chlorinity) {
        return 0.00017492711f * (0.03f + 1.805f * chlorinity);
    }

    Float downwelling_transmittance_fast(const Float &theta_i) {
        // Evaluate fitted downwelling polynomial
        // Coefficients: 0.00882913 -0.03803999  0.02232457 -0.00533469  0.97575765
        return 0.00882913 * dr::pow(theta_i, 4.0f)
                - 0.03803999 * dr::pow(theta_i, 3.0f)
                + 0.02232457 * dr::pow(theta_i, 2.0f)
                - 0.00533469 * theta_i
                + 0.97575765;
    }

    Float upwelling_transmittance_fast(const Float &theta_o) {
        // Evaluate fitted downwelling polynomial
        // Coefficients: -0.03694741  0.17051226 -0.18538372 -0.01917329  0.93060225
        return -0.03694741 * dr::pow(theta_o, 4.0f) 
            + 0.17051226 * dr::pow(theta_o, 3.0f)
            - 0.18538372 * dr::pow(theta_o, 2.0f)
            - 0.01917329 * theta_o
            + 0.93060225; 
    }

    Float downwelling_transmittance(const Float &wavelength,
                                    const Float &theta_i, const Float &phi_i,
                                    const Float &wind_direction, const Float &wind_speed, 
                                    const Float &chlorinity) {
        Float sum = 0.0f;
        Float downwelling_opacity = 0.0f;
        for (auto i = 0; i < 48; i++)
            for (auto j = 0; j < 24; j++) {
                //  Get the azimuth/zenith points and weights
                Float phi_d = azimuth_pts[i];
                Float theta_d = zenith_pts[j];
                Float azimuth_weight = azimuth_weights[i];
                Float zenith_weight = zenith_weights[j];

                // Cosine/sine of the angles
                Float c_theta_d = dr::cos(theta_d);
                Float s_theta_d = dr::sin(theta_d);

                // Compute the reflectance
                Float reflectance = eval_glint_internal(wavelength, theta_i, theta_d, phi_i, phi_d, 
                                                  wind_direction, wind_speed, chlorinity);

                // Quadrature step
                Float weight = azimuth_weight * zenith_weight;
                Float factor = c_theta_d * s_theta_d * weight;
                sum += factor;
                downwelling_opacity += reflectance * factor;
            }

        return 1.f - downwelling_opacity / sum;
    }

    Float upwelling_transmittance(const Float &wavelength,
                                  const Float &theta_o, const Float &phi_o,
                                  const Float &wind_direction, const Float &wind_speed, 
                                  const Float &chlorinity) {
        Float sum = 0.0f;
        Float upwelling_opacity = 0.0f;
        for (auto i = 0; i < 48; i++)
            for (auto j = 0; j < 24; j++) {
                //  Get the azimuth/zenith points and weights
                Float phi_d = azimuth_pts[i];
                Float theta_d = zenith_pts[j];
                Float azimuth_weight = azimuth_weights[i];
                Float zenith_weight = zenith_weights[j];

                // Cosine/sine of the angles
                Float c_theta_d = dr::cos(theta_d);
                Float s_theta_d = dr::sin(theta_d);

                // Compute the reflectance
                Float reflectance = eval_glint_internal(wavelength, theta_d, theta_o, phi_d, phi_o, 
                                                  wind_direction, wind_speed, chlorinity);

                // Quadrature step
                Float weight = azimuth_weight * zenith_weight;
                Float factor = c_theta_d * s_theta_d * weight;
                sum += factor;
                upwelling_opacity += reflectance * factor;
            }

        return 1.f - upwelling_opacity / sum;
    }

    Float downwelling_transmittance_broken(const Float &wavelength,
                                    const Float &theta_i, const Float &phi_i,
                                    const Float &wind_direction, const Float &wind_speed, 
                                    const Float &chlorinity) {
        /*
        // Gauss-Lobatto quadrature weights and nodes (preferred since bounds are included)
        std::pair<FloatStorage, FloatStorage> gs_azimuth  = quad::gauss_lobatto<FloatStorage>(48);
        std::pair<FloatStorage, FloatStorage> gs_zenith   = quad::gauss_lobatto<FloatStorage>(24);

        auto azimuth_pts = gs_azimuth.first;
        auto azimuth_weights = gs_azimuth.second;
        auto zenith_pts = gs_zenith.first;
        auto zenith_weights = gs_zenith.second;

        // Transformation of the zenith and azimuthal angles
        auto transformed_azimuth_pts = (azimuth_pts + 1.0f) * m_pi;
        auto transformed_zenith_pts = (zenith_pts + 1.0f) * (m_pi_half / 2.0f);

        // (Co)sine of the angles
        auto c_theta = dr::cos(transformed_zenith_pts);
        auto s_theta = dr::sin(transformed_zenith_pts);

        // Compute reflectance by gather ops
        for (auto i = 0; i < 48; i++)
            for (auto j = 0; j < 24; j++) {
                UInt32 azimuth_index = i;
                UInt32 zenith_index = j;

                auto test = dr::gather<Float>(transformed_azimuth_pts, azimuth_index);

                Log(Warn, "Test: %f", test);
            }
        */

        /*
        // Compute reflectance
        Float reflectance = eval_glint_internal(wavelength, theta_i, transformed_zenith_pts, phi_i, transformed_azimuth_pts,
                                       wind_direction, wind_speed, chlorinity);

        // Quadrature step
        Float weights = azimuth_weights * zenith_weights;
        Float factor = c_theta * s_theta * weights;
        Float sum = dr::sum(factor);
        Float downwelling_opacity = reflectance * factor;

        return 1.f - downwelling_opacity / sum;
        */

        return 1.0f;
    }

    Float upwelling_transmittance_broken(const Float &wavelength,
                                  const Float &theta_o, const Float &phi_o,
                                  const Float &wind_direction, const Float &wind_speed, 
                                  const Float &chlorinity) {
        /*
        // Gauss-Lobatto quadrature weights and nodes (preferred since bounds are included)
        std::pair<FloatX, FloatX> gs_azimuth  = quad::gauss_lobatto<Float>(48);
        std::pair<Float, Float> gs_zenith   = quad::gauss_lobatto<Float>(24);

        Float azimuth_pts = gs_azimuth.first;
        Float azimuth_weights = gs_azimuth.second;
        Float zenith_pts = gs_zenith.first;
        Float zenith_weights = gs_zenith.second;

        // Transformation of the zenith and azimuthal angles
        Float transformed_azimuth_pts = (azimuth_pts + 1.0f) * m_pi;
        Float transformed_zenith_pts = (zenith_pts + 1.0f) * (m_pi_half / 2.0f);
    
        // Quadrature
        Float upwelling_opacity = 0.0f;
        Float summ = 0.0f;
        for (int i = 0; i < azimuth_pts; ++i) {
            for (int j = 0; j < zenith_pts; ++j) {
                Float phi_d = transformed_azimuth_pts[i];
                Float theta_d = transformed_zenith_pts[j];
                Float azimuth_weight = azimuth_weights[i];
                Float zenith_weight = zenith_weights[j];

                // Cosine/sine of the angles
                Float c_theta_d = dr::cos(theta_d);
                Float s_theta_d = dr::sin(theta_d);

                // Compute the reflectance
                Float reflectance = eval_glint(wavelength, theta_d, theta_o, phi_d, phi_o, 
                                                  wind_direction, wind_speed, chlorinity);

                // Quadrature step
                Float weight = azimuth_weight * zenith_weight;
                Float factor = c_theta_d * s_theta_d * weight;
                summ += factor;
                upwelling_opacity += reflectance * factor;
            }
        }

        return 1.f - upwelling_opacity / summ;
        */

       return 1.0f;
    }

    Float r_omega(const Float &wavelength,
                  const Float &pigmentation) {
        Float wavelength_nm = wavelength * 1000.0f;

        // Backscattering coefficient
        Float molecular_scatter_coeff = m_ocean_props.molecular_scatter_coeff(wavelength_nm);
        Float scattering_coeff = 0.30f * dr::pow(pigmentation, 0.62);
        Float backscatter_ratio = 0.002f + 0.02f * (0.5f - 0.25f * dr::log(pigmentation)) * (550.0 / wavelength_nm);
        Float backscatter_coeff = 0.5f * molecular_scatter_coeff + scattering_coeff * backscatter_ratio;

        // (Diffuse) attenuation coefficient
        Float k = m_ocean_props.attn_k(wavelength_nm);
        Float chi = m_ocean_props.attn_chi(wavelength_nm);
        Float e = m_ocean_props.attn_e(wavelength_nm);
        Float attn_coeff = k + chi * dr::pow(pigmentation, e);

        // Iterative computation of the reflectance
        Float u = 0.75f;
        Float r_omega = 0.33f * backscatter_coeff / (u * attn_coeff);

        //  TODO: Change
        bool converged = false;
        while (!converged) {
            // Update u
            u = (0.9f * (1.0f - r_omega)) / (1.0f + 2.25f * r_omega);

            // Update reflectance
            Float r_omega_new = 0.33f * backscatter_coeff / (u * attn_coeff);

            // Create a mask that marks the converged values
            auto convergence_mask = Mask(dr::abs((r_omega_new - r_omega) / r_omega_new) < 0.001f);
            converged = dr::all(convergence_mask);

            //if (dr::abs((r_omega_new - r_omega) / r_omega_new) < 0.001f)
            //    converged = true;

            // Update reflectance ONLY for non-converged values
            r_omega = dr::select(convergence_mask, r_omega, r_omega_new);
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
        m_wavelength = props.get<ScalarFloat>("wavelength");
        m_wind_speed = props.get<ScalarFloat>("wind_speed");
        m_wind_direction = props.get<ScalarFloat>("wind_direction");
        m_chlorinity = props.get<ScalarFloat>("chlorinity");
        m_pigmentation = props.get<ScalarFloat>("pigmentation");

        // Initialize the ocean utilities
        m_ocean_utils = new OceanUtilities<Float, Spectrum>();

        // Set the BSDF flags
        // => Whitecap reflectance is "diffuse"
        m_components.push_back(BSDFFlags::DiffuseReflection | 
                               BSDFFlags::FrontSide);
    
        // => Sun glint reflectance at the water surface is "specular"
        m_components.push_back(BSDFFlags::GlossyReflection | 
                               BSDFFlags::FrontSide | BSDFFlags::BackSide);

        // => Underlight reflectance is "diffuse" but transmissive
        m_components.push_back(BSDFFlags::DiffuseTransmission | 
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

    Float eval_ocean(const Vector3f &wi, const Vector3f &wo) {
        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float whitecap_reflectance = eval_whitecaps();
        Float glint_reflectance = eval_glint(wi, wo);
        Float underlight_reflectance = eval_underlight(wi, wo);

        return (coverage * whitecap_reflectance) 
            + (1 - coverage) * glint_reflectance
            + (1 - (coverage * whitecap_reflectance)) * underlight_reflectance;
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        bool has_underlight = ctx.is_enabled(BSDFFlags::DiffuseTransmission, 2);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return { dr::zeros<BSDFSample3f>(), UnpolarizedSpectrum(0.f) };

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Vector3f wo = warp::square_to_cosine_hemisphere(sample2);

        UnpolarizedSpectrum result(0.f);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        Float coverage = m_ocean_utils->eval_whitecap_coverage(m_wind_speed);
        Float glint = std::rand() / (Float) RAND_MAX;

        // Based on the coverage, we choose to sample either the whitecaps or the sun glint / underlight
        Float prob_whitecap = coverage,
              prob_glint = (1 - coverage) * glint,
              prob_underlight = (1 - coverage) * (1.0f - glint);

        // TODO: What if one of lobes disabled?

        Mask sample_whitecap = active && (sample1 < prob_whitecap),
             sample_glint = active && (sample1 >= prob_whitecap && sample1 < prob_whitecap + prob_glint),
             sample_underlight = active && (sample1 >= prob_whitecap + prob_glint);

        if (dr::any_or<true>(sample_whitecap)) {
            dr::masked(bs.wo, sample_whitecap) = wo;
            dr::masked(bs.pdf, sample_whitecap) = warp::square_to_cosine_hemisphere_pdf(wo);
            dr::masked(bs.sampled_component, sample_whitecap) = 0;
            dr::masked(bs.sampled_type, sample_whitecap) = +BSDFFlags::DiffuseReflection; 
        }

        if (dr::any_or<true>(sample_glint)) {
            dr::masked(bs.wo, sample_glint) = wo;
            dr::masked(bs.pdf, sample_glint) = warp::square_to_cosine_hemisphere_pdf(wo);
            dr::masked(bs.sampled_component, sample_glint) = 1;
            dr::masked(bs.sampled_type, sample_glint) = +BSDFFlags::GlossyReflection;
        }

        if (dr::any_or<true>(sample_underlight)) {
            dr::masked(bs.wo, sample_underlight) = wo;
            dr::masked(bs.pdf, sample_underlight) = warp::square_to_cosine_hemisphere_pdf(wo);
            dr::masked(bs.sampled_component, sample_underlight) = 2;
            dr::masked(bs.sampled_type, sample_underlight) = +BSDFFlags::DiffuseTransmission;
        }

        bs.eta = 1.f;
        active &= bs.pdf > 0.f;
        result = eval_ocean(si.wi, wo);

        return { bs, depolarizer<Spectrum>(result) & active };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        //Log(Warn, "Evaluating OceanicBSDF");

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        bool has_underlight = ctx.is_enabled(BSDFFlags::DiffuseTransmission, 2);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);
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

        result[is_reflect] = (coverage * whitecap_reflectance) 
            + (1 - coverage) * glint_reflectance
            + (1 - (coverage * whitecap_reflectance)) * underlight_reflectance;
        
        return depolarizer<Spectrum>(result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        return dr::select(active, warp::square_to_cosine_hemisphere_pdf(wo), 0.f);
    }

    void traverse(TraversalCallback *callback) override {
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanicLegacy[" << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << std::endl
            << "  wind_direction = " << string::indent(m_wind_direction) << std::endl
            << "  chlorinity = " << string::indent(m_chlorinity) << std::endl
            << "  pigmentation = " << string::indent(m_pigmentation) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    // User-provided fields
    ScalarFloat m_wavelength;
    ScalarFloat m_wind_speed;
    ScalarFloat m_wind_direction;
    ScalarFloat m_chlorinity;
    ScalarFloat m_pigmentation;

    // Fields used to compute whitecap reflectance
    OceanUtilities<Float, Spectrum> *m_ocean_utils;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)