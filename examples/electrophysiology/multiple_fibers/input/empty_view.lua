mmCreateView("insitu", "View3D", "::testview")
mmCreateModule("SimpleSphereRenderer", "::rnd")
--mmCreateModule("adiosDataSource" "::dat")
mmCreateModule("ADIOStoMultiParticle", "::converter")

mmCreateCall("CallRender3D", "::testview::rendering", "::rnd::rendering")
