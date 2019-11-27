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
                print(pattern)
                return ["NA" for i in range(expectedNumber)]
        else:
            return info
    except:
        return ["NA" for i in range(expectedNumber)]


files = [f for f in listdir(inputDir) if isfile(join(inputDir, f))]
files = [inputDir + f for f in files if ('GSE' in f)]
print(len(files))

gsm_df = pd.DataFrame(columns=['gsm',
                               'gse',
                               'gsm_submission_date',
                               'gsm_last_update_date',
                               'organism',
                               'description',
                               'gsm_title',
                               'gsm_type',
                               'gsm_source_name',
                               'gsm_extract_protocol',
                               'gpl',
                               'gsm_library_selection',
                               'gsm_library_source',
                               'gsm_library_strategy',
                               'row_count'])

gse_df = pd.DataFrame(columns=['gse',
                               'gse_title',
                               'number_gsm',
                               'gse_gpl',
                               'gse_submission_date',
                               "gse_status",
                               'gse_last_update_date',
                               'gse_summary',
                               'gse_type',
                               'sub_series',
                               'super_series',
                               'has_single_cell',
                               'has_multichannel',
                               'mixed',
                               'gsm'])

gse_df.to_csv(gseOutFile, sep="\t", index=False)
gsm_df.to_csv(gsmOutFile, sep="\t", index=False)

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

        gse_title = re.findall(r'!Series_title\t"(.*)"', temp)[0]

        gse_status = re.findall(r'!Series_status\t"(.*)"', temp)[0]

        gse_submission_date = re.findall(r'!Series_submission_date\t"(.*)"', temp)[0]

        gse_last_update_date = re.findall(r'!Series_last_update_date\t"(.*)"', temp)[0]

        gse_summary = re.findall(r'!Series_summary\t"(.*)"', temp)[0]

        gse_type = re.findall('!Series_type\t"(.*)"', temp)
        gse_type = list(filter(None, gse_type))
        # Check if multiple types of experiments are in same matrix
        if len(list(set(gse_type))) != 1:
            mixedFlag = True
        else:
            mixedFlag = False
        gse_type = ";".join(gse_type)

        try:
            sub_series = re.findall(r'"SubSeries of: (.*)"', temp)
            sub_series = ",".join(sub_series)
        except IndexError:
            sub_series = "NA"

        try:
            super_series = re.findall(r'"SuperSeries of: (.*)"', temp)
            super_series = ",".join(super_series)
        except IndexError:
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

        data = {'gsm': gsm,
                'gse': [gse for i in range(number_gsm)],
                'gsm_submission_date': gsm_submission_date,
                'gsm_last_update_date': gsm_last_update_date,
                'organism': organism,
                'description': description,
                'gsm_title': gsm_title,
                'gsm_type': gsm_type,
                'gsm_source_name': gsm_source_name,
                'gsm_extract_protocol': gsm_extract_protocol,
                'gpl': gpl,
                'gsm_library_selection': gsm_library_selection,
                'gsm_library_source': gsm_library_source,
                'gsm_library_strategy': gsm_library_strategy,
                'row_count': row_count
                }

        df = pd.DataFrame(data)
        df.to_csv(gsmOutFile, sep="\t", mode='a', header=False, index=False)

        data = {'gse': [gse],
                'gse_title': [gse_title],
                'number_gsm': [number_gsm],
                'gse_gpl': [gse_gpl],
                'gse_submission_date': [gse_submission_date],
                "gse_status": [gse_status],
                'gse_last_update_date': [gse_last_update_date],
                'gse_summary': [gse_summary],
                'gse_type': [gse_type],
                'sub_series': [sub_series],
                'super_series': [super_series],
                'has_single_cell': [SCFlag],
                'has_multichannel': [MCFlag],
                'mixed': [mixedFlag],
                'gsm': [';'.join(gsm)]
                }
        df = pd.DataFrame(data)
        df.to_csv(gseOutFile, sep="\t", mode='a', header=False, index=False)
