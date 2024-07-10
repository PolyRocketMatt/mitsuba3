#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class OceanicBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    OceanicBSDF(const Properties &props) : Base(props) {
        // Retrieve the parameters used in 6SV
        m_wavelength = props.texture<Texture>("wavelength");
        m_wind_speed = props.texture<Texture>("wind_speed");

        // Set the BSDF flags
        // => Whitecap reflectance is "diffuse"
        m_components.push_back(BSDFFlags::DiffuseReflection | 
                               BSDFFlags::FrontSide | BSDFFlags::BackSide);
    
        // => Sun glint reflectance at the water surface

        m_flags = m_components[0];
        dr::set_attr(this, "flags", m_flags);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        /*
        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return { dr::zeros<BSDFSample3f>(), UnpolarizedSpectrum(0.f) };

        // Compute the cosine of the outgoing directons compared 
        // to the normal at the surface
        Float cos_theta_i = Frame3f::cos_theta(si.wi);

        // Sample the hemisphere with a cosine distribution
        Vector3f wo = warp::square_to_cosine_hemisphere(sample2);

        // Initialize the BSDF sample to return
        BSDFSample3f bs = zero<BSDFSample3f>();

        // Create an unpolarized spectrum to return
        UnpolarizedSpectrum value(0.f);

        // Select the lobe to be sampled
        Float whitecap_sampling_weight = 1.f;

        // Retrieve wind speeds
        UnpolarizedSpectrum wind_speed = 38.0f * m_wind_speed->eval(si, active);

        // Evaluate by activating lanes
        value = dr::select(active, Float(1.f), 0.f);

        // Koepke 1984 with linear interpolation to obtain the efficiency factor
        // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
        value[active] = 4.0f + [4.0f * (wind_speed / 38.0f) - 2.0f];

        // Compute PDF (= cosine weighed for now)
        bs.pdf = dr::select(active, warp::square_to_cosine_hemisphere_pdf(wo), 0.f);
    
        // Set other interaction fields
        bs.eta = 1.f;
        bs.sampled_component = dr::select(active, UInt32(0), UInt32(0));
        bs.sampled_type = dr::select(active, UInt32(+BSDFFlags::DiffuseReflection), 
                                             UInt32(+BSDFFlags::DiffuseReflection));

        // Return the result
        return { bs, value };
        */

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

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);

        if (has_whitecap) {
            UnpolarizedSpectrum wind_speed = m_wind_speed->eval(si, active);
            
            // Koepke 1984 with linear interpolation to obtain the efficiency factor
            // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
            auto wind_speed_frac = wind_speed / m_wind_speed_max;
            auto f_eff = m_f_eff_base + (m_f_eff_base * wind_speed_frac - m_f_eff_mod);

            // Here we compute the effective reflectance by looking up the tabulated
            // values provided by Whitlock et al. 1982
            auto index_f = (m_wavelength->eval(si, active).x() - Float(0.2)) / Float(0.1);
            auto index_i = dr::floor2int<UInt32>(index_f);
            auto eff_reflectance = dr::gather<Float>(m_eff_reflectance, index_i);

            // Based on the power-law provided by Monahan & Muircheartaigh 1980,
            // the fractional whitecap coverage can be determined
            auto coverage = m_monahan_alpha * dr::pow(wind_speed, m_monahan_lambda);

            result[active] = coverage * f_eff * eff_reflectance;
        }

        return depolarizer<Spectrum>(result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        /*
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return 0.f;

        // Ensure that incoming direction is in upper hemisphere
        Vector3f wo_flip{ wo.x(), wo.y(), dr::abs(cos_theta_o) };

        Float result = dr::select(
            active, warp::square_to_cosine_hemisphere_pdf(wo_flip), 0.f);

        return result;
        */
        return dr::select(active, warp::square_to_cosine_hemisphere_pdf(wo), 0.f);
    }

    void traverse(TraversalCallback *callback) override {
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanicLegacy[" << std::endl
            << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    // User-provided parameters
    ref<Texture> m_wavelength;
    ref<Texture> m_wind_speed;

    // Parameters used to compute whitecap reflectance
    ScalarFloat m_f_eff_base = 0.4f;
    ScalarFloat m_f_eff_mod = 0.2f;
    ScalarFloat m_wind_speed_max = 38.0f;
    ScalarFloat m_monahan_alpha = 2.951 * 10e-6;
    ScalarFloat m_monahan_lambda = 3.52;

    // These values range from 0.2 to 4.0 μm
    ScalarFloat m_eff_reflectance[39] = {
        0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210, 0.200, 0.190,
        0.175, 0.155, 0.130, 0.080, 0.100, 0.105, 0.100, 0.080, 0.045, 0.055,
        0.065, 0.060, 0.055, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
    };
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)