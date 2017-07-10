/*
A KBase module: ZS_MutualInfo
*/

module ZS_MutualInfo {
    /*
        A string representing a compound id.
    */
    typedef string compound_id;
    /* 
        The workspace ID for a FBAModel data object.
        @id ws KBaseFBA.FBAModel
    */
    typedef string ws_fbamodel_id;
    /* 
        The workspace ID for a Report object
        @id ws KBaseReport.Report
    */
	typedef string ws_report_id;
    
    typedef structure {
        ws_fbamodel_id fbamodel_id;
        list<compound_id> compounds;
    } RunFluxMutualInformationAnalysisParams;
    
    typedef structure {
        string report_name;
		ws_report_id report_ref;
    } RunFluxMutualInformationAnalysisResults;
    
    funcdef run_flux_mutual_information_analysis(RunFluxMutualInformationAnalysisParams params) returns (RunFluxMutualInformationAnalysisResults output) authentication required;
};
