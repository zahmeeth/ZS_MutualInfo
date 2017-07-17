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
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    def _generate_mutual_info(self, compounds_file, fba_file):

        df1 = pd.read_csv(fba_file)
        df1.as_matrix()


        df2 = pd.read_csv(compounds_file)
        df2.as_matrix()

        s_df1 = df1.shape
        s_df2 = df2.shape


        Temp_df2 = np.array(df2.values)
        # Create matrix with only the elements remove first column and all the rows
        Temp_df2 = Temp_df2[0:,1:]

        Bin_val_check =  np.array_equal(Temp_df2, Temp_df2.astype(bool))
        num_compounds = (s_df2[1])-1

        if ((s_df1[1]-1) != s_df2[0]) or (Bin_val_check != True) or (int(math.log(s_df2[0],2)) != num_compounds):
            print ('invalid input values')

        #-----All possible combination of the chemical compounds----------------------
        # 2.0 Sperating m0 from rest of the lables

        Temp1_df2 = df2

        cols = Temp1_df2.columns
        for i in range(1,len(cols)):
            Temp1_df2.loc[Temp1_df2[cols[i]] == 1 , cols[i]] = cols[i]

        print Temp1_df2

        # 2.1 Creating a disctionary for all FBAs except m0
        print len(Temp1_df2)
        mydict = {}
        for x in range(0,len(Temp1_df2)):
            for i in range(1,s_df2[1]):
                currentvalue = Temp1_df2.iloc[x,i]
                currentid = Temp1_df2.iloc[x,0]
                currentvalue = Temp1_df2.iloc[x,i]
                mydict.setdefault(currentid,[])
                if currentvalue > 0:
                    mydict[currentid].append(currentvalue)
                                
        # Add the first key as m0
        media_0_name = 'm0'
        mydict[media_0_name] = "['0']"
        #Sort the keys 
        mydict = collections.OrderedDict(natsort.natsorted(mydict.items())) 
        print mydict

        for k,v in mydict.iteritems():
            print k,v

        # List of Compounds combination in the list
        my_combi_list = []
        Compounds_Combi = list(range(1,num_compounds+1))
        for L in range(0, len(Compounds_Combi)+1):
          for subset in itertools.combinations(Compounds_Combi, L):
            my_combi_list.append(list(subset))
        print my_combi_list


        # Created a dictionary where the keys: 
        # list of compounds combination 
        # values are corresponding FBAs list in df2
        result_dict = {}
        for element in my_combi_list[1:]:
            for k, v in mydict.iteritems():
                if set(v).issubset(set(map(lambda x:str(x), element))):
                    key = ','.join(map(lambda x:str(x), element))
                    if result_dict.get(key):
                        media_list = result_dict[key]
                        media_list.append(k)
                        media_list = list(set(media_list))
                        result_dict.update({key: media_list})
                    else:
                        result_dict.update({key: [media_0_name, k]})            
        print result_dict

        # Created a dictionary where the keys are: 
        # list of compounds combination 
        # values are compounds combination FBAs with df1 vaules 
        All_Comp_Combi_dic = {}
        for column, value in result_dict.items():
            All_Comp_Combi_dic.update({column : df1.get(value)})   


        #To print an item from the All_Comp_Combi_dic
        df = (pd.DataFrame(All_Comp_Combi_dic.items()))

        print df[0]
        #print df[1][7]

        MI_dict = {}
        for k in range(0, len(df[0])):          
            drop_rows_df = df[1][k].drop_duplicates(keep="first")
            drop_columns_df = drop_rows_df.T.drop_duplicates(keep="first").T
            remove = []
            removed = {}
            cols = df[1][k].columns
            for i in range(len(cols)-1):
                duplicated = []
                v = df[1][k][cols[i]].values
                for j in range(i+1,len(cols)):
                    if np.array_equal(v,df[1][k][cols[j]].values):
                        remove.append(cols[j])
                        duplicated.append(cols[j])
                if duplicated and cols[i] not in remove:
                    removed.update({cols[i]:duplicated})
                count = {}
                for key, value in removed.items():
                    count.update({key: len(value)})
            
                #print v
            
                # print drop_columns_df
                values = count.values()
                # print values
                values = map(lambda x: x+1, values)
                # print values
                d = {x:values.count(x) for x in values}     

            #-------Mutual Inforamtion (MI) calculation-------------
            FBAs = len(df[1][k].columns)
            pure_entropy = math.log(FBAs,2)
            #print pure_entropy


            # If No duplicates exist and list "value" is empty
            if not values:
                #print("List is empty")
                No_duplicate_FBAs = len(drop_columns_df.columns)
                conditional_entropy = -1 * (No_duplicate_FBAs*((1/No_duplicate_FBAs)*((1/1)*math.log(1.0/1.0,2))));
                Mutual_Info = pure_entropy - conditional_entropy
                #print('Mutaul Info:', Mutual_Info)

            if values:
            # If duplicates exist and list "value" is not empty
                conditional_entropy = 0
                for key in d:
                    #print key, d[key]
                    Temp = -1 * d[key] * (key/float(FBAs)) * key * (1.0/key) * math.log(1.0/key,2)
                    conditional_entropy = Temp + conditional_entropy
                #print "%3f" %Temp
                Mutual_Info = pure_entropy - conditional_entropy
                #print('Mutaul Info:', Mutual_Info)

            MI_dict.update({df[0][k] : Mutual_Info})

            return MI_dict

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
