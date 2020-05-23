mmCreateView("adiosmegamol", "View3D", "::v")
mmCreateModule("SimpleSphereRenderer", "::rnd")
mmCreateModule("adiosDataSource","::dat")
mmCreateModule("ADIOStoMultiParticle", "::converter")
mmCreateModule("LinearTransferFunction", "::lt")

--mmCreateModule("OverrideParticleGlobals", "::override")
--mmSetParamValue("::override::overrideRadius", "True")
--mmSetParamValue("::override::radius", "0.1")

mmCreateCall("CallADIOSData", "::converter::adiosSlot", "::dat::getdata")
mmCreateCall("CallRender3D", "::v::rendering", "::rnd::rendering")
mmCreateCall("MultiParticleDataCall", "::rnd::getData", "::converter::mpSlot")
mmCreateCall("CallGetTransferFunction", "::rnd::gettransferfunction", "::lt::gettransferfunction")
