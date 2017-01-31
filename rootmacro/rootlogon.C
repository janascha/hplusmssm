void rootlogon()
{
  // Load ATLAS style
  gROOT->LoadMacro("atlas/AtlasStyle.C");
  gROOT->LoadMacro("atlas/AtlasUtils.C");
  SetAtlasStyle();
  
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
 
}
