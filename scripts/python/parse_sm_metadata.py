import re
import gzip
import sys
import pandas as pd
from os import listdir
from os.path import isfile, join

inputDir = sys.argv[1]
gseOutFile = sys.argv[2]
gsmOutFile = sys.argv[3]

print("Input directory: " + inputDir)
print("GSE out file: " + gseOutFile)
print("GSM out file: " + gsmOutFile)


# set expectedNumber = 0 if parse GSM
def parse_field(pattern, string, expectedNumber):
    try:
        info = re.findall(pattern, string)
        if len(info) > 1:
            info = [i.split("\t") for i in info]
            for i in range(1, len(info)):
                info[0] = [m + ';' + n for m, n in zip(info[0], info[i])]
            info = info[0]
        else:
            info = info[0].split("\t")
        info = [x.replace("\"", "") for x in info]
        info = list(filter(None, info))

        if expectedNumber > 0:
            if len(info) == expectedNumber:
                return info
            else:
                return ["NA" for i in range(expectedNumber)]
        else:
            return info
    except:
        return ["NA" for i in range(expectedNumber)]

gsm_df = pd.DataFrame(columns=['GSM',
                               'GSE',
                               'SUBMISSION_DATE',
                               'LAST_UPDATE_DATE',
                               'ORGANISM',
                               'DESCRIPTION',
                               'TITLE',
                               'TYPE',
                               'SOURCE_NAME',
                               'EXTRACT_PROTOCOL',
                               'GPL',
                               'LIBRARY_SELECTION',
                               'LIBRARY_SOURCE',
                               'LIBRARY_STRATEGY',
                               'ROW_COUNT'])

gse_df = pd.DataFrame(columns=['GSE',
                               'TITLE',
                               'NUMBER_GSM',
                               'GPL',
                               'SUBMISSION_DATE',
                               "STATUS",
                               'LAST_UPDATE_DATE',
                               'SUMMARY',
                               'TYPE',
                               'SUB_SERIES_OF',
                               'SUPER_SERIES_OF',
                               'HAS_SINGLE_CELL',
                               'HAS_MULTICHANNEL',
                               'MIXED'])

gse_df.to_csv(gseOutFile, sep="\t", index=False)
gsm_df.to_csv(gsmOutFile, sep="\t", index=False)

files = [f for f in listdir(inputDir) if isfile(join(inputDir, f))]
files = [inputDir + "/" + f for f in files if ('GSE' in f)]

for file in files:
    gse = re.findall('(GSE\d+\-GPL\d+|GSE\d+)', file)[0]
    print(gse)
    with gzip.open(file) as fp:
        line = fp.readline().decode("utf-8", "replace")
        temp = ""
        f = False
        expression_table = False
        count = 0
        while line:
            if ("!series_matrix_table_begin" in line):
                f = True
            if f:
                count += 1
            if count > 5:
                expression_table = True
                break
            temp += line

            line = fp.readline().decode("utf-8", "replace")

        # possible markers of SC experiment
        pattern = r"Single-cell|scRNA-seq|single cell"
        SCFlag = any(re.findall(pattern, temp, re.IGNORECASE))

        # check if experiment multichannel
        if "_ch2" in temp:
            MCFlag = True
        else:
            MCFlag = False

        # GSE
        try:
            gse_title = re.findall(r'!Series_title\t"(.*)"', temp)[0]
        except IndexError:
            next

        try:
            gse_status = re.findall(r'!Series_status\t"(.*)"', temp)[0]
        except IndexError:
            gse_status = "NA"

        try:
            gse_submission_date = re.findall(r'!Series_submission_date\t"(.*)"', temp)[0]
        except IndexError:
            gse_submission_date = "NA"

        try:
            gse_last_update_date = re.findall(r'!Series_last_update_date\t"(.*)"', temp)[0]
        except IndexError:
            gse_last_update_date = "NA"

        try:
            gse_summary = re.findall(r'!Series_summary\t"(.*)"', temp)[0]
        except IndexError:
            gse_summary = "NA"

        gse_type = re.findall('!Series_type\t"(.*)"', temp)
        gse_type = list(filter(None, gse_type))
        if len(gse_type) == 0:
            gse_type = "NA"

        # Check if multiple types of experiments are in same matrix
        if gse_type == "NA":
            mixedFlag = "NA"
        elif len(list(set(gse_type))) != 1:
            mixedFlag = True
        else:
            mixedFlag = False
        gse_type = ";".join(gse_type)

        sub_series = re.findall(r'"SubSeries of: (.*)"', temp)
        sub_series = ",".join(sub_series)
        if len(sub_series) == 0:
            sub_series = "NA"

        super_series = re.findall(r'"SuperSeries of: (.*)"', temp)
        super_series = ",".join(super_series)
        if len(super_series) == 0:
            super_series = "NA"

        # GSM
        pattern = r'!Sample_geo_accession\t(.*)'
        gsm = parse_field(pattern, temp, 0)
        number_gsm = len(gsm)

        pattern = '!Sample_title\t(.*)'
        gsm_title = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_type\t(.*)'
        gsm_type = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_submission_date\t(.*)'
        gsm_submission_date = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_last_update_date\t(.*)'
        gsm_last_update_date = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_source_name_ch1\t(.*)'
        gsm_source_name = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_molecule_ch1\t(.*)'
        gsm_molecule = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_extract_protocol_ch1\t(.*)'
        gsm_extract_protocol = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_platform_id\t(.*)'
        gpl = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_organism_ch1\t(.*)'
        organism = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_data_row_count\t(.*)'
        row_count = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_library_selection\t(.*)'
        gsm_library_selection = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_description\t(.*)'
        description = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_library_source\t(.*)'
        gsm_library_source = parse_field(pattern, temp, number_gsm)

        pattern = '!Sample_library_strategy\t(.*)'
        gsm_library_strategy = parse_field(pattern, temp, number_gsm)

        # GSE
        pattern = '!Sample_platform_id\t(.*)'
        gse_gpl = ";".join(list(set(parse_field(pattern, temp, number_gsm))))

        data = {'GSM': gsm,
                'GSE': [gse for i in range(number_gsm)],
                'SUBMISSION_DATE': gsm_submission_date,
                'LAST_UPDATE_DATE': gsm_last_update_date,
                'ORGANISM': organism,
                'DESCRIPTION': description,
                'TITLE': gsm_title,
                'TYPE': gsm_type,
                'SOURCE_NAME': gsm_source_name,
                'EXTRACT_PROTOCOL': gsm_extract_protocol,
                'GPL': gpl,
                'LIBRARY_SELECTION': gsm_library_selection,
                'LIBRARY_SOURCE': gsm_library_source,
                'LIBRARY_STRATEGY': gsm_library_strategy,
                'ROW_COUNT': row_count
                }

        df = pd.DataFrame(data)
        df.to_csv(gsmOutFile, sep="\t", mode='a', header=False, index=False)

        data = {'GSE': [gse],
                'TITLE': [gse_title],
                'NUMBER_GSM': [number_gsm],
                'GPL': [gse_gpl],
                'SUBMISSION_DATE': [gse_submission_date],
                "STATUS": [gse_status],
                'LAST_UPDATE_DATE': [gse_last_update_date],
                'SUMMARY': [gse_summary],
                'TYPE': [gse_type],
                'SUB_SERIES_OF': [sub_series],
                'SUPER_SERIES_OF': [super_series],
                'HAS_SINGLE_CELL': [SCFlag],
                'HAS_MULTICHANNEL': [MCFlag],
                'MIXED': [mixedFlag]}

        df = pd.DataFrame(data)
        df.to_csv(gseOutFile, sep="\t", mode='a', header=False, index=False)
