import csv
import plotly.express as px

def plot_d_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z=[]
        d=[]
        ang=[]
        sang=[]
        lab=[]
        ar_res=[]
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0],row[1],row[5],row[6],row[9],row[11],row[7],row[8],row[12])
            d.append(float(row[12]))
            z.append(float(row[8]))
            ang.append(float(row[17]))
            sang.append(float(row[22]))
            ar_res.append(row[11])
            lab.append(info)
    fig = px.scatter(x=d,y=z,color=ar_res,hover_name=lab,
                     labels=dict(x='Distance (Å)', y='Z-score'),)
    fig.show()


def plot_d_vs_solidangle(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z=[]
        d=[]
        ang=[]
        sang=[]
        lab=[]
        ar_res=[]
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0],row[1],row[5],row[6],row[9],row[11],row[7],row[8],row[12])
            d.append(float(row[12]))
            z.append(float(row[8]))
            ang.append(float(row[17]))
            sang.append(float(row[22]))
            ar_res.append(row[11])
            lab.append(info)
    fig = px.scatter(x=d,y=sang,color=ar_res,hover_name=lab,
                     labels=dict(x='Distance (Å)', y='Solid angle'),)
    fig.show()

def plot_azimuthal_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[5], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            d.append(float(row[12]))
            z.append(float(row[8]))
            ang.append(float(row[17]))
            sang.append(float(row[22]))
            ar_res.append(row[11])
            lab.append(info)
    fig = px.scatter(x=ang, y=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Azimuthal angle', y='Z-score'), )
    fig.show()

def plot_solidangle_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[5], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            d.append(float(row[12]))
            z.append(float(row[8]))
            ang.append(float(row[17]))
            sang.append(float(row[22]))
            ar_res.append(row[11])
            lab.append(info)
    fig = px.scatter(xs=ang, y=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Solid angle', y='Z-score'), )
    fig.show()

def plot_3d(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[5], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            d.append(float(row[12]))
            z.append(float(row[8]))
            ang.append(float(row[17]))
            sang.append(float(row[22]))
            ar_res.append(row[11])
            lab.append(info)
    fig = px.scatter_3d(x=d, y=ang, z=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Distance (Å)', z='Z-score', y='Azimuthal angle'), )
    fig.show()

