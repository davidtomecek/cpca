import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio
from pathlib import Path
import yaml

site_name = 'ikem'
run = 'resting'

match site_name:
    case 'ikem':
        config_file = Path(f"./data/{site_name}/config.yaml")
    case 'nudz':
        config_file = Path(f"./data/{site_name}/config.yaml")

with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

run_variants = config['runs'].get(run, [])

root_dir = Path('/hydra/hydra_io/vypocty/tomecek/cpca/data').joinpath(site_name)

# Load information about measured series
df_visit = pd.read_csv(next(root_dir.joinpath(run).glob('export-visit*csv')), sep=';')

column_names = df_visit.columns.to_list()

column_names = {'osobní kód (Hydra ID)': 'id', 'osobní kód (Hydra ID).1': 'hydra_id', 'Pořadí vizity': 'visit', 'Seznam sérií': 'series', 'Místo měření': 'site'}
df_visit = df_visit.rename(columns=column_names)

list_series = []
for visit in df_visit.iterrows():
    data = visit[1]

    id = data.id
    hydra_id = data.hydra_id
    visit = data.visit
    site = data.site
    series = data.series.split()

    #run_name = [s.replace(',', '') for s in series if run in s.lower()]
    run_name = [s.replace(',', '') for s in series if any(r in s.lower() for r in run_variants)]

    list_series.append({'id': id, 'visit': visit, 'run_name': run_name[0], 'site': site})

df_long = pd.DataFrame(list_series)
df_wide = df_long.pivot(index='id', columns='visit', values=['run_name', 'site']).to_dict(orient='index')

# Load information about subjects
df_person = pd.read_csv(next(root_dir.glob('export-person*csv')), skiprows=1, sep=";")

subs = df_person.iloc[:,[0,1]]
subs.columns = ['id', 'hydra_id']
visits = np.arange(1, 4)

motion = []
hydra_id = []
series = []
site = []
timepoints = []
for sub in subs.iterrows():

    id = sub[1]['id']

    fd_visits = []
    len_run_visits = []
    series_visits = []
    site_visits = []
    for visit in visits:

        try:
            sub_dir = next(root_dir.joinpath(run).glob(f"ESO*{id:05}*{visit}"))
            fd_file = next(sub_dir.glob(f"{sub_dir.name}_{run}_fd_power2012.txt"))
            fd = pd.read_csv(fd_file)
            len_run = np.size(fd, axis=0)+1
            mean_fd = fd.mean().values
            fd_visits.append(mean_fd[0])
            len_run_visits.append(len_run)
            series_visits.append(df_wide[id]['run_name', visit])
            site_visits.append(df_wide[id]['site', visit])
        except:
            fd_visits.append(np.nan)
            len_run_visits.append(np.nan)
            series_visits.append(np.nan)
            site_visits.append(np.nan)

    hydra_id.append(sub[1]['hydra_id'])
    series.append(series_visits)
    motion.append(fd_visits)
    timepoints.append(len_run_visits)
    site.append(site_visits)

df_fd = pd.concat((pd.DataFrame(hydra_id), pd.DataFrame(motion), pd.DataFrame(timepoints), pd.DataFrame(series), pd.DataFrame(site)), axis=1)
df_fd.columns = ['hydraID', 'meanFD_V1', 'meanFD_V2', 'meanFD_V3', 'length_V1', 'length_V2', 'length_V3', 'series_V1', 'series_V2', 'series_V3', 'site_V1', 'site_V2', 'site_V3']

df_fd.to_csv(root_dir.joinpath(run, f"eso_{run}_fd_motion_{site_name}.csv"), index=False)
