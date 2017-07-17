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
<<<<<<< HEAD
=======
        The workspace ID for a Media data object.
        @id ws KBaseBiochem.Media
    */
    typedef string ws_media_id;
    /* 
>>>>>>> 926dddc52326374882c0ae646eed77f45218d0f7
        The workspace ID for a Report object
        @id ws KBaseReport.Report
    */
	typedef string ws_report_id;
    
    typedef structure {
<<<<<<< HEAD
        ws_fbamodel_id fbamodel_id;
        list<compound_id> compounds;
        string workspace;
        string media_id;
=======
	ws_media_id media_id;
        ws_fbamodel_id fbamodel_id;
        list<compound_id> compounds;
>>>>>>> 926dddc52326374882c0ae646eed77f45218d0f7
    } RunFluxMutualInformationAnalysisParams;
    
    typedef structure {
        string report_name;
		ws_report_id report_ref;
    } RunFluxMutualInformationAnalysisResults;
    
    funcdef run_flux_mutual_information_analysis(RunFluxMutualInformationAnalysisParams params) returns (RunFluxMutualInformationAnalysisResults output) authentication required;
};
