#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import re


parser = argparse.ArgumentParser()
parser.add_argument('folder', nargs='?', help='folder in which pdf should be summarized')
parser.add_argument('sep', nargs='?', help='deliminator')
args = parser.parse_args()

print(args.folder)
print(args.sep)

file_names = [x for x in os.listdir(args.folder) if x.find('.pdf')]
file_information = [x.split(args.sep) for x in file_names]
file_information = pd.DataFrame(file_information)
file_information['link'] = [args.folder + '/' + x for x in file_names]
file_information = file_information.loc[:, ['link'] + list(range(file_information.shape[1] - 1))]

writer = pd.ExcelWriter(args.folder + '.xlsx', engine='xlsxwriter')
file_information.to_excel(writer, sheet_name='link')

worksheet = writer.sheets['link']
for i in range(file_information.shape[0]):
    x = re.sub('/var/www/html/', 'http://35.237.17.250/', file_information['link'][i])
    worksheet.write_url('B' + str(i + 2), x)

writer.save()
writer.close()
