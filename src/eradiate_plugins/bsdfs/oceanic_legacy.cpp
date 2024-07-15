#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

// Header content - Potentially move somewhere else later
#ifndef OCEAN_PROPS
#define OCEAN_PROPS

template<typename Float, typename Spectrum>
class OceanProperties {
public:
    MI_IMPORT_TYPES()

    OceanProperties(const Properties &props) {
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

        m_ior_real = IrregularContinuousDistribution<Float>(
            ior_wavelengths.data(), ior_real_data.data(), ior_real_data.size()
        );

        m_ior_imag = IrregularContinuousDistribution<Float>(
            ior_wavelengths.data(), ior_cplx_data.data(), ior_cplx_data.size()
        );
    }

private:
    //  2. Complex IOR of water (Hale & Querry 1973)
    IrregularContinuousDistribution<Float> m_ior_real;
    IrregularContinuousDistribution<Float> m_ior_imag;
};

template<typename Float, typename Spectrum>
class CoxMunkDistribution {
public:
    MI_IMPORT_TYPES()
private:
};

template<typename Float, typename Spectrum>
class OceanUtilities {
public:
    MI_IMPORT_TYPES()

    Spectrum eval_cox_munk(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi_i, const Spectrum &phi_w, const Spectrum &wind_speed) const {
        // Difference between the azimuth of th incoming light direction and 
        // the azimuth of the wind speed.
        Spectrum chi = phi_i - phi_w;

        // Compute the Cox-Munk slope distribution, term by term
        Spectrum s_c = dr::sqrt(cm_s_c_sqr(wind_speed));
        Spectrum s_u = dr::sqrt(cm_s_u_sqr(wind_speed));

        Spectrum ksi = cm_z_x_prime(theta_i, theta_o, chi) / s_c;
        Spectrum eta = cm_z_y_prime(theta_i, theta_o, chi) / s_u;
    
        Spectrum ksi_sqr = ksi * ksi;
        Spectrum eta_sqr = eta * eta;

        Spectrum normalization_c = 1.f / (dr::TwoPi<Float> * s_c * s_u);
        Spectrum exp_factor = dr::exp(-0.5 * (ksi_sqr + eta_sqr));

        Spectrum a = (cm_c_21(wind_speed) / 2.0f) * (ksi_sqr - 1.0f) * eta;
        Spectrum b = (cm_c_03(wind_speed) / 6.0f) * (eta_sqr * eta - 3.0f * eta);
        Spectrum c = (m_c_40 / 24.0f) * (ksi_sqr * ksi_sqr - 6.0f * ksi_sqr + 3.0f);
        Spectrum d = (m_c_22 / 4.0f) * (ksi_sqr - 1.0f) * (eta_sqr - 1.0f);
        Spectrum e = (m_c_04 / 24.0f) * (eta_sqr * eta_sqr - 6.0f * eta_sqr + 3.0f);

        // Combine
        return normalization_c * exp_factor * (1.0f - a - b + c + d + e);
    }

    Spectrum eval_fresnel(const Spectrum &theta_i, const Spectrum &theta_o, 
                  const Spectrum &phi_i, const Spectrum &phi_o, 
                  const Spectrum &n_real, const Spectrum &n_cplx,
                  const Spectrum &chlorinity) {
        Spectrum phi_rel = phi_i - phi_o;
        Spectrum chi = fresnel_chi(theta_i, theta_o, phi_rel);
        Spectrum cos_chi = dr::cos(chi);

        //  Compute the real and complex parts of the index of refraction 
        //  with Friedman and Sverdrup correction
        Spectrum salinity = friedman_sverdrup_salinity(chlorinity);
        Spectrum n_real_corrected = n_real + m_salinity_factor * salinity;
        Spectrum n_real_sqr = dr::sqr(n_real_corrected);
        Spectrum n_cplx_sqr = dr::sqr(n_cplx);

        //  Compute a values
        Spectrum a1 = fresnel_a_1(n_real_sqr, n_cplx_sqr, chi);
        Spectrum a2 = fresnel_a_2(n_real_sqr, n_cplx_sqr, chi);

        //  Compute b values
        Spectrum b1 = vb_1(n_real_sqr, n_cplx_sqr, chi);
        Spectrum b2 = fresnel_b_2(n_real_corrected, n_cplx, chi);

        //  Compute u² and v²
        Spectrum u_sqr = fresnel_u_sqr(a1, a2);
        Spectrum v_sqr = fresnel_v_sqr(a1, a2);

        //  Compute u and v
        Spectrum u = dr::sqrt(u_sqr);
        Spectrum v = dr::sqrt(v_sqr);

        //  Compute the Fresnel reflection coefficient
        Spectrum left = (dr::sqr(b1 - u) + dr::sqr(b2 + v)) / (dr::sqr(b1 + u) + dr::sqr(b2 - v));
        Spectrum right = (dr::sqr(cos_chi - u) + v_sqr) / (dr::sqr(cos_chi + u) + v_sqr);

        return 0.5f * (left + right);
    }
private:
    // Cox-Munk distribution parameters
    ScalarFloat m_c_40 = 0.40f;
    ScalarFloat m_c_22 = 0.12f;
    ScalarFloat m_c_04 = 0.23f;

    // Fresnel parameters
    ScalarFloat m_salinity_factor = 0.00017492711f;

    // Underlight parameters

    Spectrum cm_c_21(const Spectrum &wind_speed) const {
        return 0.01f - 0.0086f * wind_speed;
    }

    Spectrum cm_c_03(const Spectrum &wind_speed) const {
        return 0.04f - 0.033f * wind_speed;
    }

    Spectrum cm_s_c_sqr(const Spectrum &wind_speed) const {
        return 0.003f + 0.00192f * wind_speed;
    }

    Spectrum cm_s_u_sqr(const Spectrum &wind_speed) const {
        return 0.00316f * wind_speed;
    }

    Spectrum cm_z_x(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi) const {
        return (-dr::sin(theta_o) * dr::sin(phi)) / (dr::cos(theta_i) * dr::cos(theta_o));
    }

    Spectrum cm_z_y(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi) const {
        return (dr::sin(theta_i) + dr::sin(theta_o) * dr::cos(phi)) / (dr::cos(theta_i) * dr::cos(theta_o));
    }

    Spectrum cm_z_x_prime(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &chi) const {
        return cm_z_x(theta_i, theta_o, chi) * dr::cos(chi) + cm_z_y(theta_i, theta_o, chi) * dr::sin(chi);
    }

    Spectrum cm_z_y_prime(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &chi) const {
        return -cm_z_x(theta_i, theta_o, chi) * dr::sin(chi) + cm_z_y(theta_i, theta_o, chi) * dr::cos(chi);
    }

    // Correction to the IOR of water according to Friedman (1969) and Sverdrup (1942)
    Spectrum friedman_sverdrup_salinity(const Spectrum &chlorinity) {
        return 0.03f + 1.805f * chlorinity;
    }

    Spectrum fresnel_chi(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi) const {
        auto cos_two_chi = dr::cos(theta_o) * dr::cos(theta_i) + dr::sin(theta_o) * dr::sin(theta_i) * dr::cos(phi);
        return dr::acos(cos_two_chi) / 2.0f;
    }

    Spectrum fresnel_u_sqr(const Spectrum &a_1, const Spectrum &a_2) {
        return dr::abs(a_1 + a_2) / 2.0f;
    }

    Spectrum fresnel_v_sqr(const Spectrum &a_1, const Spectrum &a_2) {
        return dr::abs(-a_1 + a_2) / 2.0f;
    }

    Spectrum fresnel_a_1(const Spectrum &n_real_sqr, const Spectrum &n_cplx_sqr, const Spectrum &chi) {
        return dr::abs(n_real_sqr - n_cplx_sqr - dr::sqr(dr::sin(chi)));
    }

    Spectrum fresnel_a_2(const Spectrum &n_real_sqr, const Spectrum &n_cplx_sqr, const Spectrum &chi) {
        Spectrum t_1 = (n_real_sqr - n_cplx_sqr - dr::sqr(dr::sin(chi)));
        Spectrum t_2 = 4.0f * n_real_sqr * n_cplx_sqr;
        return dr::sqrt(dr::sqr(t_1) + t_2);
    }

    Spectrum fresnel_b_1(const Spectrum &n_real_sqr, const Spectrum &n_cplx_sqr, const Spectrum &chi) {
        return (n_real_sqr - n_cplx_sqr) * dr::cos(chi);
    }

    Spectrum fresnel_b_2(const Spectrum &n_real, const Spectrum &n_cplx, const Spectrum &chi) {
        return 2.0f * (n_real + n_cplx) * dr::cos(chi);
    }

    Spectrum underlight_downwelling_transmittance() {
        // Integrate both azimuthal and zenithal angles
        // TODO: Make this a Monte Carlo estimtor
        Spectrum opacity = 0.0f;

        for (UInt32 zenith = 0; zenith < 90; zenith += 1) {
                // To radian
                Float theta = dr::deg_to_rad(zenith);

                // Construct a vector for the azimuthal angle. 
                // To maximize JIT, the azimuthal angle can be
                // computed in parallel.
                Float phis = dr::linspace(Float(0), Float(360), 360);
                Float cos_theta = dr::cos(theta);
                Float sin_theta = dr::sin(theta);

            }

        return 1.0f - opacity;
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
        m_wavelength = props.texture<Texture>("wavelength");
        m_wind_speed = props.texture<Texture>("wind_speed");
        m_wind_direction = props.texture<Texture>("wind_direction");
        m_salinity = props.texture<Texture>("salinity");

        //  Wavelengths considered by 6SV
        std::vector<ScalarFloat> wc_wavelengths = { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
                                                    2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1,
                                                    3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0 };

        // Effective reflectance values for whitecaps, given originally
        // by Whitlock et al. 1982
        std::vector<ScalarFloat> wc_data =    { 0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210, 0.200, 0.190,
                                                0.175, 0.155, 0.130, 0.080, 0.100, 0.105, 0.100, 0.080, 0.045, 0.055,
                                                0.065, 0.060, 0.055, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };

        m_eff_reflectance = IrregularContinuousDistribution<Float>(
            wc_wavelengths.data(), wc_data.data(), wc_data.size()
        );

        // Set the BSDF flags
        // => Whitecap reflectance is "diffuse"
        m_components.push_back(BSDFFlags::DiffuseReflection | 
                               BSDFFlags::FrontSide | BSDFFlags::BackSide);
    
        // => Sun glint reflectance at the water surface
        m_components.push_back(BSDFFlags::GlossyReflection | 
                               BSDFFlags::FrontSide | BSDFFlags::BackSide);

        m_flags = m_components[0];
        dr::set_attr(this, "flags", m_flags);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        // Test
        UnpolarizedSpectrum value(0.f);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        bs.pdf = dr::select(active, warp::square_to_cosine_hemisphere_pdf(si.wi), 0.f);
        bs.eta = 1.f;
        bs.sampled_component = dr::select(active, UInt32(0), UInt32(0));
        bs.sampled_type = dr::select(active, UInt32(+BSDFFlags::DiffuseReflection), 
                                             UInt32(+BSDFFlags::DiffuseReflection));
    
        return { bs, value };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Log(Warn, "Evaluating OceanicBSDF");

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);
        
        // Compute the wind speed
        Spectrum wind_speed = m_wind_speed->eval(si, active) * m_max_wind_speed;

        // Based on the power-law provided by Monahan & Muircheartaigh 1980,
        // the fractional whitecap coverage of whitecaps can be determined.
        // Only if whitecaps are present, we compute the coverage.
        Spectrum coverage = has_whitecap ? m_monahan_alpha * dr::pow(wind_speed, m_monahan_lambda) : 0.0f;

        if (has_whitecap) {
            // Koepke 1984 with linear interpolation to obtain the efficiency factor
            // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
            Spectrum efficiency = m_f_eff_base;

            // Here we compute the effective reflectance by looking up the tabulated
            // values provided by Whitlock et al. 1982
            Float eff_reflectance = m_eff_reflectance.eval_pdf(m_wavelength->eval(si, active).x(), active);      

            // Compute the whitecap reflectance
            Spectrum whitecap_reflectance = coverage * efficiency * eff_reflectance;
        
            // TODO: Multiply by the cosine?
            // Add the whitecap reflectance to the result
            result[active] = whitecap_reflectance;
        } 

        if (has_glint) {
            // To compute the solar glint, we need the incoming and outgoing
            // zenithal angles, the incoming azimuthal angle and the wind direction
            // which is given by the wind direction texture.
            Spectrum wind_direction = m_wind_direction->eval(si, active);
            Float sin_theta_i = Frame3f::sin_theta(si.wi),
                  sin_theta_o = Frame3f::sin_theta(wo),
                  cos_phi_i = Frame3f::cos_phi(si.wi);

            // Transform sines and cosine into angles
            Float theta_i = dr::asin(sin_theta_i),
                  theta_o = dr::asin(cos_theta_o),
                  phi_i = dr::acos(cos_phi_i);

            // Using the provided data, we can compute the probability of having
            // a solar specular reflection using the Cox-Munk distribution.
            Spectrum specular_prob = m_cox_munk.eval(theta_i, theta_o, phi_i, wind_direction, wind_speed);
        
            Log(Warn, "Specular probability: %s", specular_prob);
        }

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
            << "  m_wind_direction = " << string::indent(m_wind_direction) << std::endl
            << "  salinity = " << string::indent(m_salinity) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    // User-provided fields
    ref<Texture> m_wavelength;
    ref<Texture> m_wind_speed;
    ref<Texture> m_wind_direction;
    ref<Texture> m_salinity;

    // Fields used to compute whitecap reflectance
    ScalarFloat m_f_eff_base = 0.4f;
    ScalarFloat m_monahan_alpha = 2.951f * 1e-6;
    ScalarFloat m_monahan_lambda = 3.52f;
    ScalarFloat m_max_wind_speed = 37.241869f;

    // Fields used to compute the sun glint reflectance
    CoxMunkDistribution<Float, Spectrum> m_cox_munk;

    // Distributions used for
    //  1. Effective reflectance of whitecaps (Whitlock et al 1982, Koepke 1984)
    IrregularContinuousDistribution<Float> m_eff_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)