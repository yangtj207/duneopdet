
cet_build_plugin(PDSPCalPhoton art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardataobj::RecoBase
  FFTW3::FFTW3
)

cet_build_plugin(PDSMeanRMSDT art::tool
  LIBRARIES PRIVATE
  ROOT::Hist
  #ROOT_BASIC_LIB_LIST
)

cet_build_plugin(GausOpHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::HitFinder
  larreco::CandidateHitFinderTool
  larreco::PeakFitterTool
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Utilities
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::MathCore
  TBB::tbb
)

install_headers()
install_fhicl()
install_source()
