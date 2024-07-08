#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class Oceanic final : public BSDF<Float, Spectrum> {
public:
private:
};

MI_IMPLEMENT_CLASS_VARIANT(Oceanic, BSDF)
MI_EXPORT_PLUGIN(Oceanic, "Oceanic material")
NAMESPACE_END(mitsuba)