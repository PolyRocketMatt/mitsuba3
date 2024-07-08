#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class Oceanic final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    Oceanic(const Properties &props) : Base(props) {
        // Retrieve the parameters used in 6SV
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
        UnpolarizedSpectrum wind_speed = m_wind_speed->eval(si, active);

        // Evaluate by activating lanes
        value = dr::select(active, Float(1.f), 0.f);

        // Koepke 1984 with linear interpolation to obtain the efficiency factor
        // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
        value[active] = 4.0 + [4.0 * (wind_speed / 38.0) - 2.0];

        // Compute PDF (= cosine weighed for now)
        bs.pdf = dr::select(active, warp::square_to_cosine_hemisphere_pdf(wo), 0.f);
    
        // Set other interaction fields
        bs.eta = 1.f;
        bs.sampled_component = dr::select(active, UInt32(0), UInt32(0));
        bs.sampled_type = dr::select(active, UInt32(+BSDFFlags::DiffuseReflection), 
                                             UInt32(+BSDFFlags::DiffuseReflection));

        // Return the result
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

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);

        if (has_whitecap) {
            UnpolarizedSpectrum wind_speed = m_wind_speed->eval(si, active);

            // Koepke 1984 with linear interpolation to obtain the efficiency factor
            // The maximum wind speed is chosen to be 38 m/s as to not exceed fractional coverage limits
            result[active] = 4.0 + [4.0 * (wind_speed / 38.0) - 2.0];
        }

        // TODO: Multiply by the cosine term???

        return dr::select(active, result, 0.f);
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_whitecap = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        if (unlikely(dr::none_or<false>(active) || !has_whitecap))
            return 0.f;

        // Ensure that incoming direction is in upper hemisphere
        Vector3f wo_flip{ wo.x(), wo.y(), dr::abs(cos_theta_o) };

        Float result = dr::select(
            active, warp::square_to_cosine_hemisphere_pdf(wo_flip), 0.f);

        return result;
    }

    void traverse(TraversalCallback *callback) override {
    }

    std::string to_string() const override {
        std::ostringstream oss;
        return oss.str();
    }
private:
    // Paremeters to compute fractional whitecap coverage based on 
    // "optimal power-law descrption of oceanic whitecap coverage 
    // dependence on wind speed" by Monahan et al.
    ref<Texture> m_wind_speed;
    float m_monahan_alpha = 2.951 * 10e-6;
    float m_monahan_lambda = 3.52
};

MI_IMPLEMENT_CLASS_VARIANT(Oceanic, BSDF)
MI_EXPORT_PLUGIN(Oceanic, "Oceanic material")
NAMESPACE_END(mitsuba)