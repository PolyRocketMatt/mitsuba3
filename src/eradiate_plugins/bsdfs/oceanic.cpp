#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class Oceanic final : public BSDF<Float, Spectrum> {
public:
    Oceanic(const Properties &props) : Base(props) {
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
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