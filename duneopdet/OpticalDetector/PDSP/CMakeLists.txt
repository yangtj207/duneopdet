
cet_build_plugin(PDSPCalPhoton art::EDProducer
  LIBRARIES PRIVATE
    lardataobj::RecoBase
    FFTW3::FFTW3
)

cet_build_plugin(PDSMeanRMSDT art::tool
  LIBRARIES PRIVATE
    ROOT_BASIC_LIB_LIST
)

install_headers()
install_fhicl()
install_source()
