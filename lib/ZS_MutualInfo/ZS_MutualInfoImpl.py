# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import pandas as pd
import numpy as np
import csv
import math
import collections
import natsort
import itertools
import operator
from itertools import islice
#END_HEADER


class ZS_MutualInfo:
    '''
    Module Name:
    ZS_MutualInfo

    Module Description:
    A KBase module: ZS_MutualInfo
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/zahmeeth/ZS_MutualInfo.git"
    GIT_COMMIT_HASH = "4d78d4c2ce5dc7e311c170babdc62c6461d56ea9"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass


    def run_flux_mutual_information_analysis(self, ctx, params):
        """
        :param params: instance of type
           "RunFluxMutualInformationAnalysisParams" -> structure: parameter
           "fbamodel_id" of type "ws_fbamodel_id" (The workspace ID for a
           FBAModel data object. @id ws KBaseFBA.FBAModel), parameter
           "compounds" of list of type "compound_id" (A string representing a
           compound id.), parameter "workspace" of String, parameter
           "media_id" of String
        :returns: instance of type "RunFluxMutualInformationAnalysisResults"
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of type "ws_report_id" (The workspace ID for a Report
           object @id ws KBaseReport.Report)
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_flux_mutual_information_analysis

        self.validate_params(params)

        fbamodel_id = params.get('fbamodel_id')
        compounds = params.get('compounds')
        media_id = params.get('media_id')
        workspace_name = params.get('workspaces')

        compounds_file = self._generate_compounds_file(comounds)
        fba_file = self._generate_fba_file(fba_object_ref)
        mutual_info = self._generate_mutual_info(compounds_file, fba_file)

        #END run_flux_mutual_information_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_flux_mutual_information_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
