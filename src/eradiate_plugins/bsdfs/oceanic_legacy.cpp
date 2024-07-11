#include <mitsuba/core/distr_1d.h>
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

        // Set the internal parameters on initialization
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
            UnpolarizedSpectrum wind_speed = m_wind_speed->eval(si, active) * m_max_wind_speed;

            // Koepke 1984 with linear interpolation to obtain the efficiency factor
            // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
            Spectrum efficiency = 0.4; //m_f_eff_base + (m_f_eff_base * wind_speed - m_f_eff_mod);

            // Here we compute the effective reflectance by looking up the tabulated
            // values provided by Whitlock et al. 1982
            Float eff_reflectance = m_eff_reflectance.eval_pdf(m_wavelength->eval(si, active).x(), active);      

            // Based on the power-law provided by Monahan & Muircheartaigh 1980,
            // the fractional whitecap coverage can be determined
            Spectrum coverage = m_monahan_alpha * dr::pow(wind_speed, m_monahan_lambda);

            Spectrum ref_whitecap = coverage * efficiency * eff_reflectance;

            result[active] = ref_whitecap;// * cos_theta_o;
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
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    // User-provided parameters
    ref<Texture> m_wavelength;
    ref<Texture> m_wind_speed;
    //Bool m_efficiency_correction;

    // Parameters used to compute whitecap reflectance
    ScalarFloat m_f_eff_base = 0.4f;
    ScalarFloat m_f_eff_mod = 0.2f;
    ScalarFloat m_monahan_alpha = 2.951 * 1e-6;
    ScalarFloat m_monahan_lambda = 3.52;
    ScalarFloat m_max_wind_speed = 37.241869f;

    // These values range from 0.2 to 4.0 μm
    IrregularContinuousDistribution<Float> m_eff_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanicBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanicBSDF, "Oceanic material")
NAMESPACE_END(mitsuba)