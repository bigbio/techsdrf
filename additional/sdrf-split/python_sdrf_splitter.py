import csv

import pandas as pd
import click

@click.command()
@click.option('--input', help='sdrf')
@click.option('--property', help='property to split')
@click.option('--subproperty', help='Use a subproperty for spliting')
@click.option('--project', help= 'project accession to be added to the sdrf file name')
def main(input, property, subproperty, project):
    """Simple program that greets NAME for a total of COUNT times."""
    df = pd.read_csv(input, sep='\t', skip_blank_lines=False)
    name_project = ""
    if(project is not None):
      name_project = project + "-"
    if(subproperty is None):
      d = dict(tuple(df.groupby(property)))
    else:
      d = dict(tuple(df.groupby([property, subproperty])))
    for key in d:
      dataframe = d[key]
      if(subproperty is None):
        name_file = name_project + key.replace(" ", "-") + '.sdrf.tsv'
      else:
        name_file = name_project + key[0].replace(" ", "-") + "-" + key[1].replace(" ", "-") + '.sdrf.tsv'
      dataframe.to_csv(name_file ,sep='\t', quoting=csv.QUOTE_NONE, index=False)
    print(d)

if __name__ == '__main__':
    main()
