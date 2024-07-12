#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

// Header content - Potentially move somewhere else later
#ifndef COX_MUNK_H
#define COX_MUNK_H

template<typename Float, typename Spectrum>
class CoxMunkDistribution {
public:
    MI_IMPORT_TYPES()

    Spectrum eval(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi_i, const Spectrum &phi_w, const Spectrum &wind_speed) const {
        // Difference between the azimuth of th incoming light direction and 
        // the azimuth of the wind speed.
        Spectrum chi = phi_i - phi_w;

        // Compute the Cox-Munk slope distribution, term by term
        Spectrum s_c = dr::sqrt(s_c_sqr(wind_speed));
        Spectrum s_u = dr::sqrt(s_u_sqr(wind_speed));

        Spectrum ksi = z_x_prime(theta_i, theta_o, chi) / s_c;
        Spectrum eta = z_y_prime(theta_i, theta_o, chi) / s_u;
    
        Spectrum ksi_sqr = ksi * ksi;
        Spectrum eta_sqr = eta * eta;

        Spectrum normalization_c = 1.f / (dr::TwoPi<Float> * s_c * s_u);
        Spectrum exp_factor = dr::exp(-0.5 * (ksi_sqr + eta_sqr));

        Spectrum a = (c_21(wind_speed) / 2.0f) * (ksi_sqr - 1.0f) * eta;
        Spectrum b = (c_03(wind_speed) / 6.0f) * (eta_sqr * eta - 3.0f * eta);
        Spectrum c = (m_c_40 / 24.0f) * (ksi_sqr * ksi_sqr - 6.0f * ksi_sqr + 3.0f);
        Spectrum d = (m_c_22 / 4.0f) * (ksi_sqr - 1.0f) * (eta_sqr - 1.0f);
        Spectrum e = (m_c_04 / 24.0f) * (eta_sqr * eta_sqr - 6.0f * eta_sqr + 3.0f);

        // Combine
        return normalization_c * exp_factor * (1.0f - a - b + c + d + e);
    }
private:
    ScalarFloat m_c_40 = 0.40f;
    ScalarFloat m_c_22 = 0.12f;
    ScalarFloat m_c_04 = 0.23f;

    Spectrum c_21(const Spectrum &wind_speed) const {
        return 0.01f - 0.0086f * wind_speed;
    }

    Spectrum c_03(const Spectrum &wind_speed) const {
        return 0.04f - 0.033f * wind_speed;
    }

    Spectrum s_c_sqr(const Spectrum &wind_speed) const {
        return 0.003f + 0.00192f * wind_speed;
    }

    Spectrum s_u_sqr(const Spectrum &wind_speed) const {
        return 0.00316f * wind_speed;
    }

    Spectrum z_x(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi) const {
        return (-dr::sin(theta_o) * dr::sin(phi)) / (dr::cos(theta_i) * dr::cos(theta_o));
    }

    Spectrum z_y(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &phi) const {
        return (dr::sin(theta_i) + dr::sin(theta_o) * dr::cos(phi)) / (dr::cos(theta_i) * dr::cos(theta_o));
    }

    Spectrum z_x_prime(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &chi) const {
        return z_x(theta_i, theta_o, chi) * dr::cos(chi) + z_y(theta_i, theta_o, chi) * dr::sin(chi);
    }

    Spectrum z_y_prime(const Spectrum &theta_i, const Spectrum &theta_o, const Spectrum &chi) const {
        return -z_x(theta_i, theta_o, chi) * dr::sin(chi) + z_y(theta_i, theta_o, chi) * dr::cos(chi);
    }
};

#endif // COX_MUNK_H

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

        // Effective reflectance values for whitecaps, given originally
        // by Whitlock et al. 1982
        std::vector<ScalarFloat> wc_wavelengths = { 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f,
                                                    1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.1f,
                                                    2.2f, 2.3f, 2.4f, 2.5f, 2.6f, 2.7f, 2.8f, 2.9f, 3.0f, 3.1f,
                                                    3.2f, 3.3f, 3.4f, 3.5f, 3.6f, 3.7f, 3.8f, 3.9f, 4.0f };
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

    // These values range from 0.2 to 4.0 μm
    IrregularContinuousDistribution<Float> m_eff_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)