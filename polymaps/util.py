import re, io
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import re

mass2element = {
    16: 'O',
    12: "C",
    1: "H",
    14: 'N', 
    -999: 'X',
}
element2color = {
    "O": 'rgba(255, 0, 0, 1.0)',
    "C": 'rgba(105, 105, 105, 1.0)',
    "H": 'rgba(255, 255, 255, 1.0)',
    'N': 'rgba(255, 0, 255, 1.0)',
    "X": 'rgba(138, 43, 226, 1.0)',
}
element2size = {
    "O": 74.,
    'N': 74.,
    "C": 77.,
    "H": 46.,
    "X": 90.,
}



def read_lmp_template(fname):
    lmp_str = ""
    with open(fname, 'r') as lmp_data:
        lmp_str = lmp_data.read()

    sec_names = []
    for sec_name in re.finditer('\n+[A-Z][a-z]*.*\n+', lmp_str):
        sec_names.append(sec_name.group(0))

    sec_strs = []
    next_str = lmp_str

    for sec_name in sec_names:
        _str_list = next_str.split(sec_name)
        sec_strs.append(_str_list[0])
        next_str = _str_list[1]


    sec_names = ['head',] + sec_names
    sec_strs.append(next_str)
    
    return sec_names, sec_strs

def write_lmp_template(fname, sec_names, sec_strs):
    
    lmp_str = sec_strs[0]
    
    for i in range(1, len(sec_strs)):
        lmp_str = lmp_str + '\n\n' + sec_names[i] + '\n\n' + sec_strs[i]
    
    with open(fname, 'w') as lmp_data:
        lmp_data.write(lmp_str)
    
    print('Written to file: ' + fname)
    return

def viz_lmp_template(filename):
    sec_names, sec_strs =read_lmp_template(filename)
    double_bonds=[]
    triple_bonds=[]
    sec_df = {}
    
    for i in range(1, len(sec_names)):
        
        if 'Types' in sec_names[i]:
            names = ['id', 'type']
        elif 'Charges' in sec_names[i]:
            names = ['id', 'q']
        elif 'Coords' in sec_names[i]:
            names = ['id', 'x', 'y', 'z']
        elif 'Bonds' in sec_names[i]:
            names = ['id', 'type', 'at1', 'at2']
        elif 'Angles' in sec_names[i]:
            names = ['id', 'type', 'at1', 'at2', 'at3']
        elif 'Dihedrals' in sec_names[i]:
            names = ['id', 'type', 'at1', 'at2', 'at3', 'at4']
        elif 'Impropers' in sec_names[i]:
            names = ['id', 'type', 'at1', 'at2', 'at3', 'at4']

        df = pd.read_csv(io.StringIO(sec_strs[i]), sep=r'\s+', names=names).reset_index(drop=True)
        df.index = df.index + 1
        sec_df[sec_names[i]] = df.copy(deep=True)
        del df
        df=None
    df = pd.concat([sec_df['\n\nTypes\n\n'], sec_df['\n\nCharges\n\n'], sec_df['\n\nCoords\n\n'], ], axis=1)
    df['mol'] = 1
    df = df[['id', 'mol', 'type', 'q', 'x', 'y', 'z']]
    df = df.loc[:,~df.columns.duplicated()]

    df['color'] = 'rgba(155, 155, 0, 1.0)'
    df['size'] = 15.2*3
    df['tip']=None
    bond_df = sec_df['\n\nBonds\n\n']
    bond_df['bond_order'] = 1
    bond_df.loc[double_bonds, 'bond_order'] = 2
    bond_df.loc[triple_bonds, 'bond_order'] = 3
    viz_mol(df, bond_df, annotation=True)

def read_lmp_data(fname):
    lmp_str = ""
    with open(fname, 'r') as lmp_data:
        lmp_str = lmp_data.read()
    
    sec_names = []
    for sec_name in re.finditer('\n+[A-Z][a-z]*.*\n+', lmp_str):
        sec_names.append(sec_name.group(0))
    
    sec_strs = []
    next_str = lmp_str
    mass_sec_id = -1
    atom_sec_id = -1
    bond_sec_id = -1
    angle_sec_id = -1
    dihedral_sec_id = -1
    improper_sec_id = -1
    
    lj_coeff_sec_id = -1
    bond_coeff_sec_id = -1
    angle_coeff_sec_id = -1
    dihedral_coeff_sec_id = -1
    improper_coeff_sec_id = -1
    for sec_name in sec_names:
        _str_list = next_str.split(sec_name)
        sec_strs.append(_str_list[0])
        next_str = _str_list[1]
        if "Masses" in sec_name:
            mass_sec_id = sec_names.index(sec_name) + 1
        if "Atoms" in sec_name:
            atom_sec_id = sec_names.index(sec_name) + 1
        if "Bonds" in sec_name:
            bond_sec_id = sec_names.index(sec_name) + 1
        if "Angles" in sec_name:
            angle_sec_id = sec_names.index(sec_name) + 1
        if "Dihedrals" in sec_name:
            dihedral_sec_id = sec_names.index(sec_name) + 1
        if "Impropers" in sec_name:
            improper_sec_id = sec_names.index(sec_name) + 1
        if 'Pair Coeffs' in sec_name:
            lj_coeff_sec_id = sec_names.index(sec_name) + 1
        if 'Bond Coeffs' in sec_name:
            bond_coeff_sec_id = sec_names.index(sec_name) + 1
        if 'Angle Coeffs' in sec_name:
            angle_coeff_sec_id = sec_names.index(sec_name) + 1
        if 'Dihedral Coeffs' in sec_name:
            dihedral_coeff_sec_id = sec_names.index(sec_name) + 1
        if 'Improper Coeffs' in sec_name:
            improper_coeff_sec_id = sec_names.index(sec_name) + 1

    sec_names = ['head',] + sec_names
    sec_strs.append(next_str)
    
    return sec_names, sec_strs, mass_sec_id, atom_sec_id, bond_sec_id, angle_sec_id, dihedral_sec_id, improper_sec_id, lj_coeff_sec_id, bond_coeff_sec_id, angle_coeff_sec_id, dihedral_coeff_sec_id, improper_coeff_sec_id

def write_lmp_data(fname, sec_names, sec_strs):
    lmp_str = sec_strs[0]
    for i in range(1, len(sec_strs)):
        lmp_str = lmp_str + sec_names[i] + sec_strs[i]
    with open(fname, 'w') as lmp_data:
        lmp_data.write(lmp_str)
    print('Written to file: ' + fname)
    return

def viz_mol(_df, _bond_df, annotation=False, box=None, writeHTML=False):
    
    df = _df.copy(deep=True)
    bond_df = _bond_df.copy(deep=True)
    
    df.loc[:, 'text'] = ['atom-'] * len(df)
    df.loc[:, 'text'] = df.loc[:, 'text'] + df.loc[:, 'id'].astype('str')

    df['tip'] = 'id: ' + df['id'].astype('str') + '<br>type: ' + df['type'].astype('str')
    if "element" in df.columns:
        df['tip'] = df['tip'] + '<br>element: ' + df['element'].astype('str')
    if "q" in df.columns:
        df['tip'] = df['tip'] + '<br>charge: ' + df['q'].astype('str')
    if "mol" in df.columns:
        df['tip'] = df['tip'] + '<br>mol: ' + df['mol'].astype('str')
    if "comment" in df.columns:
        df['tip'] = df['tip'] + '<br>comment: ' + df['comment'].astype('str')
    bond_df.loc[:, 'text'] = ['bond-'] * len(bond_df)
    bond_df.loc[:, 'text'] = bond_df.loc[:, 'text'] + bond_df.loc[:, 'id'].astype('str')
    
    bond_dict_x = df.set_index('id').to_dict()['x']
    bond_dict_y = df.set_index('id').to_dict()['y']
    bond_dict_z = df.set_index('id').to_dict()['z']
    bond_df.loc[:, 'x1'] = bond_df['at1'].map(bond_dict_x)
    bond_df.loc[:, 'y1'] = bond_df['at1'].map(bond_dict_y)
    bond_df.loc[:, 'z1'] = bond_df['at1'].map(bond_dict_z)

    bond_df.loc[:, 'x2'] = bond_df['at2'].map(bond_dict_x)
    bond_df.loc[:, 'y2'] = bond_df['at2'].map(bond_dict_y)
    bond_df.loc[:, 'z2'] = bond_df['at2'].map(bond_dict_z)
    
    df.loc[:, 'showarrow'] = False
    df.loc[:, 'opacity'] = 0.8
    df.loc[:, 'font']=[{'color':'rgba(0,0,255,1)',}] * len(df)
    bond_df.loc[:, 'showarrow'] = False
    bond_df.loc[:, 'opacity'] = 0.8
    bond_df.loc[:, 'font']=[{'color':'rgba(0,0,255,1)',}] * len(bond_df)
    bond_df.loc[:, 'x'] = (bond_df.loc[:, 'x1'] + bond_df.loc[:, 'x2']) * 0.5
    bond_df.loc[:, 'y'] = (bond_df.loc[:, 'y1'] + bond_df.loc[:, 'y2']) * 0.5
    bond_df.loc[:, 'z'] = (bond_df.loc[:, 'z1'] + bond_df.loc[:, 'z2']) * 0.5

    
    if annotation:
        annotation_list = df[['x', 'y', 'z', 'showarrow', 'opacity', 'text', 'font']].to_dict('records') + \
            bond_df[['x', 'y', 'z', 'showarrow', 'opacity', 'text', 'font']].to_dict('records')
    else:
        annotation_list = []
    data=[go.Scatter3d(x=df['x'], y=df['y'], z=df['z'], 
                                       text=df['tip'],
                                       mode='markers',
                                       marker=dict(
                                           color=df['color'],
                                           size=df['size'] * 3.4 / (len(df)**(1./3.)),
                                           opacity=1.0,
                                       ))]
    
    for index, row in bond_df.iterrows():
        data.append(go.Scatter3d(x=[row['x1'], row['x2']], 
                                 y=[row['y1'], row['y2']], 
                                 z=[row['z1'], row['z2']], 
                                 mode='lines',
                                 line=dict(
                                     color='rgba(220, 200, 200, 0.5)',
                                     width=5,
                                 )))
        if row['bond_order'] == 2:
            
            
            data.append(go.Scatter3d(x=[row['x1'], row['x2']], 
                                     y=[row['y1'], row['y2']], 
                                     z=[row['z1'], row['z2']],
                                     mode='lines',
                                     line=dict(
                                         color='rgba(150, 150, 150, 0.35)',
                                         width=30,
                                     )))
            
            if row['bond_order'] == 3:            
                data.append(go.Scatter3d(x=[row['x1'], row['x2']], 
                                         y=[row['y1'], row['y2']], 
                                         z=[row['z1'], row['z2']],
                                         mode='lines',
                                         line=dict(
                                             color='rgba(100, 100, 100, 0.25)',
                                             width=50,
                                         )))
                
    xyzmin = min([df['x'].min(), df['y'].min(), df['z'].min()])
    xyzmax = max([df['x'].max(), df['y'].max(), df['z'].max()])
    if box != None:
        xlo = box[0][0]
        xhi = box[0][1]
        ylo = box[1][0]
        yhi = box[1][1]
        zlo = box[2][0]
        zhi = box[2][1]
        
        xyzmin = min([df['x'].min(), df['y'].min(), df['z'].min(), xlo, ylo, zlo])
        xyzmax = max([df['x'].max(), df['y'].max(), df['z'].max(), xhi, yhi, zhi])
        annotation_list.append({'x': xhi, 
                                'y': ylo, 
                                'z': zlo, 
                                'showarrow': False, 
                                'opacity': 1.0, 
                                'text': "X",  
                                'font': {'size': 30, 'color':'rgba(255,0,0,1)',}})
        annotation_list.append({'x': xlo, 
                                'y': yhi, 
                                'z': zlo, 
                                'showarrow': False, 
                                'opacity': 1.0, 
                                'text': "Y",  
                                'font': {'size': 30, 'color':'rgba(0,255,0,1)',}})
        annotation_list.append({'x': xlo, 
                                'y': ylo, 
                                'z': zhi, 
                                'showarrow': False, 
                                'opacity': 1.0, 
                                'text': "Z",  
                                'font': {'size': 30, 'color':'rgba(0,0,255,1)',}})
        data.append(go.Scatter3d(x=[xlo, xhi], 
                                 y=[ylo, ylo], 
                                 z=[zlo, zlo], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(225, 0, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xlo], 
                                 y=[ylo, yhi], 
                                 z=[zlo, zlo], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(0, 225, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xlo], 
                                 y=[ylo, ylo], 
                                 z=[zlo, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(0, 0, 225, 1.0)',
                                     width=3,
                                 )))
        
        data.append(go.Scatter3d(x=[xhi, xhi], 
                                 y=[ylo, ylo], 
                                 z=[zlo, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xlo], 
                                 y=[yhi, yhi], 
                                 z=[zlo, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xhi, xhi], 
                                 y=[yhi, yhi], 
                                 z=[zlo, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        
        data.append(go.Scatter3d(x=[xlo, xhi], 
                                 y=[ylo, ylo], 
                                 z=[zhi, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xhi], 
                                 y=[yhi, yhi], 
                                 z=[zlo, zlo], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xhi], 
                                 y=[yhi, yhi], 
                                 z=[zhi, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        
        data.append(go.Scatter3d(x=[xhi, xhi], 
                                 y=[ylo, yhi], 
                                 z=[zlo, zlo], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xlo, xlo], 
                                 y=[ylo, yhi], 
                                 z=[zhi, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        data.append(go.Scatter3d(x=[xhi, xhi], 
                                 y=[ylo, yhi], 
                                 z=[zhi, zhi], 
                                 mode='lines', 
                                 line=dict(
                                     color='rgba(255, 255, 0, 1.0)',
                                     width=3,
                                 )))
        
    
    fig = go.Figure(data=data)
    
    
    DeltaX = xyzmax - xyzmin
    fig.update_layout(
        scene = dict(
            annotations=annotation_list,
            xaxis = dict(nticks=10, range=[xyzmin-10,xyzmax+10],
                         backgroundcolor="rgba(80, 70, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            yaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 80, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            zaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 70, 80, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                         ),
        ),
        width=1400,
        height=1400,
        margin=dict(r=10, l=10, b=10, t=10),
        showlegend=False)
    fig.update_layout(scene_aspectmode='cube', paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')
    fig.update_layout(
        scene_aspectmode='cube', 
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig.update_layout(
        font_color="rgba(150,150,150,1)",
        title_font_color="rgba(150,150,150,1)",
        legend_title_font_color="rgba(150,150,150,1)",
    )

    if writeHTML:
        fig.write_html("vizmol.html")
    fig.show()
    
    
    
def open_lmp_data(lammps_data_file, viz=True, annotation=False, double_bonds=[], triple_bonds=[], box=False, unwrap=False, resize=False):
    
    sec_names, sec_strs, mass_sec_id, atom_sec_id, bond_sec_id, angle_sec_id, dihedral_sec_id, improper_sec_id, lj_coeff_sec_id, bond_coeff_sec_id, angle_coeff_sec_id, dihedral_coeff_sec_id, improper_coeff_sec_id = read_lmp_data(lammps_data_file)
    
    L_df = pd.read_csv(io.StringIO('\n'.join([x for x in sec_strs[0].split('\n') if "xlo" in x or "ylo" in x or "zlo" in x])), 
                sep=r'\s+', 
                names=['_min', '_max'], 
                usecols=[0, 1]).reset_index(drop=True)

    xhi = L_df.loc[0, '_max']
    xlo = L_df.loc[0, '_min']
    yhi = L_df.loc[1, '_max']
    ylo = L_df.loc[1, '_min']
    zhi = L_df.loc[2, '_max']
    zlo = L_df.loc[2, '_min']

    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo

    mass_df = pd.read_csv(io.StringIO(sec_strs[mass_sec_id]), sep=r'\s+', names=['type', 'mass'], usecols=[0, 1]).reset_index(drop=True)
    mass_df['element'] = mass_df['mass'].round().map(mass2element)
    type_element_map = mass_df.set_index('type').to_dict()['element']
    
    
    if unwrap:
        try:
            df = pd.read_csv(io.StringIO(sec_strs[atom_sec_id]), sep=r'\s+', 
                             names=['id', 'mol', 'type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz'], 
                             usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]).reset_index(drop=True)
        except:
            df = pd.read_csv(io.StringIO(sec_strs[atom_sec_id]), sep=r'\s+', 
                             names=['id', 'mol', 'type', 'q', 'x', 'y', 'z'], 
                             usecols=[0, 1, 2, 3, 4, 5, 6]).reset_index(drop=True)
            df["ix"] = 0
            df["iy"] = 0
            df["iz"] = 0
        df.index = df.index + 1
        
        df['_x'] = df['x']        
        df['x'] = df['x'] + (Lx * df['ix'])
        df['_y'] = df['y']        
        df['y'] = df['y'] + (Ly * df['iy'])
        df['_z'] = df['z']        
        df['z'] = df['z'] + (Lx * df['iz'])
    else:
        df = pd.read_csv(io.StringIO(sec_strs[atom_sec_id]), sep=r'\s+', names=['id', 'mol', 'type', 'q', 'x', 'y', 'z'], usecols=[0, 1, 2, 3, 4, 5, 6]).reset_index(drop=True)
        df.index = df.index + 1
        if resize:
            xhi = df["x"].max() + 1.
            xlo = df["x"].min() - 1.
            yhi = df["y"].max() + 1.
            ylo = df["y"].min() - 1.
            zhi = df["z"].max() + 1.
            zlo = df["z"].min() - 1.
            Lx = xhi - xlo
            Ly = yhi - ylo
            Lz = zhi - zlo

    df['element'] = df['type'].map(type_element_map)
    df['color'] = df['element'].map(element2color)
    df['size'] = df['element'].map(element2size)
    
    try:
        bond_df = pd.read_csv(io.StringIO(sec_strs[bond_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2'], usecols=[0, 1, 2, 3]).reset_index(drop=True)
        bond_df.index = bond_df.index + 1
        bond_df['bond_order'] = 1
        bond_df.loc[double_bonds, 'bond_order'] = 2
        bond_df.loc[triple_bonds, 'bond_order'] = 3
    except:
        bond_df = pd.DataFrame()
        pass

    try:
        angle_df = pd.read_csv(io.StringIO(sec_strs[angle_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3'], usecols=[0, 1, 2, 3, 4]).reset_index(drop=True)
        angle_df.index = angle_df.index + 1
    except:
        angle_df = pd.DataFrame()
        pass
    
    
    try:
        dihedral_df = pd.read_csv(io.StringIO(sec_strs[dihedral_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3', 'at4'], usecols=[0, 1, 2, 3, 4, 5]).reset_index(drop=True)
        dihedral_df.index = dihedral_df.index + 1
    except:
        dihedral_df = pd.DataFrame()
        pass
    
    try:
        improper_df = pd.read_csv(io.StringIO(sec_strs[improper_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3', 'at4'], usecols=[0, 1, 2, 3, 4, 5]).reset_index(drop=True)
        improper_df.index = improper_df.index + 1
    except:
        improper_df = pd.DataFrame()
        pass
    
    lj_coeff = None
    if lj_coeff_sec_id != -1:
        lj_coeff = pd.read_csv(io.StringIO(sec_strs[lj_coeff_sec_id]), names=['type', 'epsilon', 'sigma'], sep=r'\s+')
    bond_coeff = None
    if bond_coeff_sec_id != -1:
        bond_coeff = pd.read_csv(io.StringIO(sec_strs[bond_coeff_sec_id]), names=['type', 'kr', 'r0'], sep=r'\s+')
    angle_coeff = None
    if angle_coeff_sec_id != -1:
        angle_coeff = pd.read_csv(io.StringIO(sec_strs[angle_coeff_sec_id]), names=['type', 'ktheta', 'theta0'], sep=r'\s+')
    dihedral_coeff = None
    if dihedral_coeff_sec_id != -1:
        dihedral_coeff = pd.read_csv(io.StringIO(sec_strs[dihedral_coeff_sec_id]), names=['type', 'k1', 'k2', 'k3', 'k4'], sep=r'\s+')
    improper_coeff = None
    if improper_coeff_sec_id != -1:
        improper_coeff = pd.read_csv(io.StringIO(sec_strs[improper_coeff_sec_id]), names=['type', 'k', 'd', 'n'], sep=r'\s+')
    if viz==True: 
        if box==True:
            viz_mol(df, bond_df, annotation=annotation, box=[[xlo, xhi], [ylo, yhi], [zlo, zhi]])
        else:
            viz_mol(df, bond_df, annotation=annotation, box=None)
    
    return mass_df, df, bond_df, angle_df, dihedral_df, improper_df, [[xlo, xhi], [ylo, yhi], [zlo, zhi]], lj_coeff, bond_coeff, angle_coeff, dihedral_coeff, improper_coeff



def read_xyz_traj(filename='traj.xyz'):
    colnames = ["at", "x", "y", "z"]
    with open(filename, 'r') as myfile:
        xyz = myfile.read()
    frame_natoms = []
    frame_nsteps = []
    frame_xyzs = []
    line_list = xyz.split('\n')
    i = 0
    while i < len(line_list):
        if line_list[i].isdigit():
            frame_natoms.append(int(line_list[i]))
            frame_nsteps.append(int(re.sub("[^0-9]", "", line_list[i+1])))
            frame_xyzs.append(pd.read_csv(io.StringIO("\n".join(line_list[i+2:i+2+frame_natoms[-1]])), sep=r"\s+", names=colnames))
            i = i + 2 + frame_natoms[-1]
        else:
            i = i + 1
    print("Read xyz file with: " + str(frame_natoms[-1]) + " atoms, " + str(len(frame_nsteps)) + " frames. ")
    return frame_natoms, frame_nsteps, frame_xyzs


def viz_xyz_traj(_natoms, _nsteps, _xyzs, mass_df, writeHTML=False):
    _xyz_lims = pd.DataFrame([[_xyzs[i]["x"].min(), _xyzs[i]["x"].max(), 
                               _xyzs[i]["y"].min(), _xyzs[i]["y"].max(), 
                               _xyzs[i]["z"].min(), _xyzs[i]["z"].max()] for i in range(0, len(_xyzs))], 
                             columns=["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"])
    xyzmin = min(_xyz_lims["xmin"].min(), _xyz_lims["ymin"].min(), _xyz_lims["zmin"].min())
    xyzmax = max(_xyz_lims["xmax"].max(), _xyz_lims["ymax"].max(), _xyz_lims["zmax"].max())
    duration = 150

    # make figure
    fig_dict = {
        "data": [],
        "layout": {},
        "frames": []
    }

    # fill in most of layout
    fig_dict["layout"]["hovermode"] = "closest"
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": duration, 
                                              "redraw": True},
                                    "fromcurrent": True, 
                                    "mode": "immediate",
                                    "transition": {"duration": duration,
                                                   "easing": "quadratic-in-out"}
                                   }],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": duration, 
                                                "redraw": False},
                                      "mode": "immediate",
                                      "transition": {"duration": duration}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Step:",
            "visible": True,
            "xanchor": "left"
        },
        "transition": {"duration": duration, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }


    # make frames
    for k in range(0, len(_xyzs)):
        _xyzs[k]["element"] = _xyzs[k]["at"].map(dict(zip(mass_df["type"], mass_df["element"])))
        _xyzs[k]["mass"] = _xyzs[k]["at"].map(dict(zip(mass_df["type"], mass_df["mass"])))
        _xyzs[k]["color"] = _xyzs[k]["element"].map(element2color)
        _xyzs[k]["size"] = _xyzs[k]["element"].map(element2size)
        _xyzs[k]['id'] = _xyzs[k].index + 1
        _xyzs[k]['tip'] = "id: " + _xyzs[k]['id'].astype('str') + '<br>' + 'type: ' + _xyzs[k]['at'].astype('str') + '<br>'
        curr_frame = go.Frame(
            data=[
                go.Scatter3d(
                #dict(
                    mode="markers",
                    x=_xyzs[k]["x"],
                    y=_xyzs[k]["y"],
                    z=_xyzs[k]["z"],
                    text=_xyzs[k]['tip'],
                    marker=dict(
                        color=_xyzs[k]["color"],
                        size=_xyzs[k]["size"] * 3.4 / (len(_xyzs[k])**(1./3.)),
                        opacity=1.0,
                    ),
                )
            ],
            name=str(_nsteps[k]),
        )


        fig_dict["frames"].append(curr_frame)

        slider_step = {
            "args": [
                [str(_nsteps[k])],
                {
                    "frame": {
                        "duration": duration,
                        "redraw": True,
                    },
                    "mode": "immediate",
                    "transition": {
                        "duration": duration,
                    }
                }
            ],
            "label": str(_nsteps[k]),
            "method": "animate",
            "name": str(_nsteps[k]),
        }
        sliders_dict["steps"].append(slider_step)


    fig_dict['data'] = fig_dict["frames"][0]['data']
    fig_dict["layout"]["sliders"] = [sliders_dict]
    fig = go.Figure(fig_dict)

    DeltaX = xyzmax - xyzmin
    annotation_list = []
    fig.update_layout(
        scene = dict(
            annotations=annotation_list,
            xaxis = dict(nticks=10, range=[xyzmin-10,xyzmax+10],
                         backgroundcolor="rgba(80, 70, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            yaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 80, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            zaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 70, 80, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                         ),
        ),
        width=1400,
        height=1400,
        margin=dict(r=10, l=10, b=10, t=10),
        showlegend=False)
    fig.update_layout(scene_aspectmode='cube', paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')
    fig.update_layout(
        scene_aspectmode='cube', 
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig.update_layout(
        font_color="rgba(150,150,150,1)",
        title_font_color="rgba(150,150,150,1)",
        legend_title_font_color="rgba(150,150,150,1)",
    )

    if writeHTML:
        fig.write_html("xyztraj.html")
    fig.show()
    
def lmp_traj2pandas(fname="traj.lammpstrj"):
    traj = None
    with open(fname, 'r') as myfile:
        traj = myfile.read()

    _nsteps = []
    _natoms = []
    _xyzs = []
    traj_list = traj.split("\n")

    xlo = []
    xhi = []
    ylo = []
    yhi = []
    zlo = []
    zhi = []
    i = 0
    while i < len(traj_list):
        if traj_list[i].startswith("ITEM:"):
            if "TIMESTEP" in traj_list[i]:
                _nsteps.append(int(traj_list[i+1]))
                _natoms.append(int(traj_list[i+3]))
                _x, _y, _z = traj_list[i+5:i+8]
                _xlo, _xhi = _x.split(" ")
                _ylo, _yhi = _y.split(" ")
                _zlo, _zhi = _z.split(" ")
                xlo.append(float(_xlo))
                xhi.append(float(_xhi))
                ylo.append(float(_ylo))
                yhi.append(float(_yhi))
                zlo.append(float(_zlo))
                zhi.append(float(_zhi))
                _xyz = pd.read_csv(io.StringIO(traj_list[i+8].replace("ITEM: ATOMS ", "") + "\n" + "\n".join(traj_list[i+9:i+9+_natoms[-1]])), header=0, sep=r"\s+")
                _xyzs.append(_xyz)
                i = i+9+_natoms[-1]
        else:
            i = i + 1
    return _xyzs

def viz_lmp_traj(type2cgenff, fname="traj.lammpstrj", writeHTML=False):
    traj = None
    with open(fname, 'r') as myfile:
        traj = myfile.read()

    _nsteps = []
    _natoms = []
    _xyzs = []
    traj_list = traj.split("\n")


    xlo = []
    xhi = []
    ylo = []
    yhi = []
    zlo = []
    zhi = []
    i = 0
    while i < len(traj_list):
        if traj_list[i].startswith("ITEM:"):
            if "TIMESTEP" in traj_list[i]:
                _nsteps.append(int(traj_list[i+1]))
                _natoms.append(int(traj_list[i+3]))
                _x, _y, _z = traj_list[i+5:i+8]
                _xlo, _xhi = _x.split(" ")
                _ylo, _yhi = _y.split(" ")
                _zlo, _zhi = _z.split(" ")
                xlo.append(float(_xlo))
                xhi.append(float(_xhi))
                ylo.append(float(_ylo))
                yhi.append(float(_yhi))
                zlo.append(float(_zlo))
                zhi.append(float(_zhi))
                _xyz = pd.read_csv(io.StringIO(traj_list[i+8].replace("ITEM: ATOMS ", "") + "\n" + "\n".join(traj_list[i+9:i+9+_natoms[-1]])), header=0, sep=r"\s+")
                _xyz["element"] = _xyz["mass"].round(0).map(mass2element)
                _xyz["cgenff"] = _xyz["type"].map(type2cgenff)
                _xyz["color"] = _xyz["element"].map(element2color)
                _xyz["size"] = _xyz["element"].map(element2size)
                _xyz['tip'] = "id: " + _xyz['id'].astype('str') + '<br>' + \
                              'type: ' + _xyz['type'].astype('str') + '<br>' + \
                              'cgenff: ' + _xyz['cgenff'].astype('str') + '<br>' + \
                              'charge: ' + _xyz["q"].apply(lambda i: ("+" if i > 0 else "") + str(i)).astype('str')
                _xyzs.append(_xyz)
                i = i+9+_natoms[-1]
        else:
            i = i + 1


    xyzmin = min(min(xlo), min(ylo), min(zlo))
    xyzmax = max(max(xhi), max(yhi), max(zhi))

    duration = 100

    # make figure
    fig_dict = {
        "data": [],
        "layout": {},
        "frames": []
    }

    # fill in most of layout
    fig_dict["layout"]["hovermode"] = "closest"
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": duration, 
                                              "redraw": True},
                                    "fromcurrent": True, 
                                    "mode": "immediate",
                                    "transition": {"duration": duration,
                                                   "easing": "quadratic-in-out"}
                                   }],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": duration, 
                                                "redraw": False},
                                      "mode": "immediate",
                                      "transition": {"duration": duration}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Step:",
            "visible": True,
            "xanchor": "left"
        },
        "transition": {"duration": duration, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }

    annotation_list = []
    annotation_list.append({'x': xhi[0]+1, 
                            'y': ylo[0], 
                            'z': zlo[0], 
                            'showarrow': False, 
                            'opacity': 1.0, 
                            'text': "X",  
                            'font': {'size': 30, 'color':'rgba(255,0,0,1)',}})
    annotation_list.append({'x': xlo[0], 
                            'y': yhi[0]+1, 
                            'z': zlo[0], 
                            'showarrow': False, 
                            'opacity': 1.0, 
                            'text': "Y",  
                            'font': {'size': 30, 'color':'rgba(0,255,0,1)',}})
    annotation_list.append({'x': xlo[0], 
                            'y': ylo[0], 
                            'z': zhi[0]+1, 
                            'showarrow': False, 
                            'opacity': 1.0, 
                            'text': "Z",  
                            'font': {'size': 30, 'color':'rgba(0,0,255,1)',}})
    # make frames
    for k in range(0, len(_xyzs)):

        _xyzs[k]["color"] = _xyzs[k]["element"].map(element2color)
        _xyzs[k]["size"] = _xyzs[k]["element"].map(element2size)
        curr_frame = go.Frame(
            data=[
                go.Scatter3d(
                    mode="markers",
                    x=_xyzs[k]["x"],
                    y=_xyzs[k]["y"],
                    z=_xyzs[k]["z"],
                    text=_xyzs[k]['tip'],
                    marker=dict(
                        color=_xyzs[k]["color"],
                        size=_xyzs[k]["size"] * 3.4 / (len(_xyzs[k])**(1./3.)),
                        opacity=1.0,
                    ),
                ),

                go.Scatter3d(x=[xlo[k], xhi[k]], 
                    y=[ylo[k], ylo[k]], 
                    z=[zlo[k], zlo[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(225, 0, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xlo[k]], 
                    y=[ylo[k], yhi[k]], 
                    z=[zlo[k], zlo[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(0, 225, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xlo[k]], 
                    y=[ylo[k], ylo[k]], 
                    z=[zlo[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(0, 0, 225, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xhi[k], xhi[k]], 
                    y=[ylo[k], ylo[k]], 
                    z=[zlo[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xlo[k]], 
                    y=[yhi[k], yhi[k]], 
                    z=[zlo[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xhi[k], xhi[k]], 
                    y=[yhi[k], yhi[k]], 
                    z=[zlo[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xhi[k]], 
                    y=[ylo[k], ylo[k]], 
                    z=[zhi[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xhi[k]], 
                    y=[yhi[k], yhi[k]], 
                    z=[zlo[k], zlo[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xhi[k]], 
                    y=[yhi[k], yhi[k]], 
                    z=[zhi[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xhi[k], xhi[k]], 
                    y=[ylo[k], yhi[k]], 
                    z=[zlo[k], zlo[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xlo[k], xlo[k]], 
                    y=[ylo[k], yhi[k]], 
                    z=[zhi[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                )), 
                go.Scatter3d(x=[xhi[k], xhi[k]], 
                    y=[ylo[k], yhi[k]], 
                    z=[zhi[k], zhi[k]], 
                    mode='lines', 
                    line=dict(
                     color='rgba(255, 255, 0, 1.0)',
                     width=3,
                ))
            ],
            name=str(_nsteps[k]),
        )






        fig_dict["frames"].append(curr_frame)

        slider_step = {
            "args": [
                [str(_nsteps[k])],
                {
                    "frame": {
                        "duration": duration,
                        "redraw": True,
                    },
                    "mode": "immediate",
                    "transition": {
                        "duration": duration,
                    }
                }
            ],
            "label": str(_nsteps[k]),
            "method": "animate",
            "name": str(_nsteps[k]),
        }
        sliders_dict["steps"].append(slider_step)


    fig_dict['data'] = fig_dict["frames"][0]['data']
    fig_dict["layout"]["sliders"] = [sliders_dict]
    fig = go.Figure(fig_dict)

    DeltaX = xyzmax - xyzmin
    fig.update_layout(
        scene = dict(
            annotations=annotation_list,
            xaxis = dict(nticks=10, range=[xyzmin-10,xyzmax+10],
                         backgroundcolor="rgba(80, 70, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            yaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 80, 70, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                        ),
            zaxis = dict(nticks=10, range=[xyzmin-10, xyzmax+10],
                         backgroundcolor="rgba(70, 70, 80, 0.5)",
                         gridcolor="white",
                         showbackground=True,
                         zerolinecolor="white",
                         title=dict(font=dict(color="rgba(150,150,150,1)")),
                         ),
        ),
        width=1400,
        height=1400,
        margin=dict(r=10, l=10, b=10, t=10),
        showlegend=False)
    fig.update_layout(scene_aspectmode='cube', paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')
    fig.update_layout(
        scene_aspectmode='cube', 
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig.update_layout(
        font_color="rgba(150,150,150,1)",
        title_font_color="rgba(150,150,150,1)",
        legend_title_font_color="rgba(150,150,150,1)",
    )
    if writeHTML:
        fig.write_html(fname+".html")
    fig.show()


def viz_lmp_data(lammps_data_file, annotation=False, double_bonds=[], triple_bonds=[]):
    #lammps_data_file = './mcmd/ec/monomer.data'
    sec_names, sec_strs, mass_sec_id, atom_sec_id, bond_sec_id, angle_sec_id, dihedral_sec_id, improper_sec_id, lj_coeff_sec_id, bond_coeff_sec_id, angle_coeff_sec_id, dihedral_coeff_sec_id, improper_coeff_sec_id = read_lmp_data(lammps_data_file)
    #sec_names, sec_strs, mass_sec_id, atom_sec_id, bond_sec_id, angle_sec_id, dihedral_sec_id, improper_sec_id = read_lmp_data(lammps_data_file)

    mass_df = pd.read_csv(io.StringIO(sec_strs[mass_sec_id]), sep=r'\s+', names=['type', 'mass'], usecols=[0, 1]).reset_index(drop=True)
    mass_df['element'] = mass_df['mass'].round().map({
        16: 'O',
        12: "C",
        1: "H",
        14: 'N', 
    })
    type_element_map = mass_df.set_index('type').to_dict()['element']

    df = pd.read_csv(io.StringIO(sec_strs[atom_sec_id]), sep=r'\s+', names=['id', 'mol', 'type', 'q', 'x', 'y', 'z'], usecols=[0, 1, 2, 3, 4, 5, 6]).reset_index(drop=True)
    df.index = df.index + 1
    df['element'] = df['type'].map(type_element_map)
    df['color'] = df['element'].map({
        "O": 'rgba(255, 255, 255, 1.0)',
        "C": 'rgba(105, 105, 105, 1.0)',
        "H": 'rgba(255, 0, 0, 1.0)',
        'N': 'rgba(255, 0, 255, 1.0)',
    })
    df['size'] = df['element'].map({
        "O": 15.2*3,
        'N': 16.3*3,
        "C": 17.0*3,
        "H": 12.0*3,
    })
    df['<br>'] = '<br>'
    df['id: '] = 'id: '
    df['type: '] = 'type: '
    df['charge: '] = 'charge: '
    df['mol: '] = 'mol: '
    df['tip'] = df['id: '] + df['id'].astype('str') + df['<br>'] + df['type: '] + df['type'].astype('str') + df['<br>'] + df['mol: '] + df['mol'].astype('str') + df['<br>'] + df['charge: '] + df['q'].astype('str')
    
    bond_df = pd.read_csv(io.StringIO(sec_strs[bond_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2'], usecols=[0, 1, 2, 3]).reset_index(drop=True)
    bond_df.index = bond_df.index + 1
    bond_df['bond_order'] = 1
    bond_df.loc[double_bonds, 'bond_order'] = 2
    bond_df.loc[triple_bonds, 'bond_order'] = 3
    #print(df, bond_df)
    angle_df = pd.read_csv(io.StringIO(sec_strs[angle_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3'], usecols=[0, 1, 2, 3, 4]).reset_index(drop=True)
    angle_df.index = angle_df.index + 1
    
    dihedral_df = pd.read_csv(io.StringIO(sec_strs[dihedral_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3', 'at4'], usecols=[0, 1, 2, 3, 4, 5]).reset_index(drop=True)
    dihedral_df.index = dihedral_df.index + 1
    
    improper_df = pd.read_csv(io.StringIO(sec_strs[improper_sec_id]), sep=r'\s+', names=['id', 'type', 'at1', 'at2', 'at3', 'at4'], usecols=[0, 1, 2, 3, 4, 5]).reset_index(drop=True)
    improper_df.index = improper_df.index + 1
    
    viz_mol(df, bond_df, annotation=annotation)
    return mass_df, df, bond_df, angle_df, dihedral_df, improper_df
    
def build_react_template(df, bond_df, angle_df, dihedral_df, improper_df, selected_atoms, template_file_name):
    atom_id_remap = dict(zip(selected_atoms, [x + 1 for x in range(len(selected_atoms))]))

    out_type = df.loc[selected_atoms, 'type'].reset_index(drop=True).copy(deep=True)
    out_type.index += 1
    out_q = df.loc[selected_atoms, 'q'].reset_index(drop=True).copy(deep=True)
    out_q.index += 1
    out_xyz = df.loc[selected_atoms, ['x', 'y', 'z']].reset_index(drop=True).copy(deep=True)
    out_xyz.index += 1
    out_str = "\n\nTypes\n\n" + out_type.to_string(header=False) + "\n\nCharges\n\n" + out_q.to_string(header=False) + "\n\nCoords\n\n" + out_xyz.to_string(header=False)
    
    bond_set1 = set(bond_df[bond_df['at1'].isin(selected_atoms)].index.to_list())
    bond_set2 = set(bond_df[bond_df['at2'].isin(selected_atoms)].index.to_list())
    bond_list = list(bond_set1 & bond_set2)
    bond_df2 = bond_df.loc[bond_list, :].reset_index(drop=True).copy(deep=True)
    bond_df2.index = bond_df2.index + 1
    
    bond_df2 = bond_df2.assign(at1=bond_df2.loc[:,'at1'].map(atom_id_remap))
    bond_df2 = bond_df2.assign(at2=bond_df2.loc[:,'at2'].map(atom_id_remap))
    out_str = out_str + "\n\nBonds\n\n" + bond_df2[['type', 'at1', 'at2']].to_string(header=False)
    
    angle_set1 = set(angle_df[angle_df['at1'].isin(selected_atoms)].index.to_list())
    angle_set2 = set(angle_df[angle_df['at2'].isin(selected_atoms)].index.to_list())
    angle_set3 = set(angle_df[angle_df['at3'].isin(selected_atoms)].index.to_list())
    angle_list = list(angle_set1 & angle_set2 & angle_set3)
    angle_df2 = angle_df.loc[angle_list, :].reset_index(drop=True).copy(deep=True)
    angle_df2.index = angle_df2.index + 1
    angle_df2 = angle_df2.assign(at1=angle_df2.loc[:,'at1'].map(atom_id_remap))
    angle_df2 = angle_df2.assign(at2=angle_df2.loc[:,'at2'].map(atom_id_remap))
    angle_df2 = angle_df2.assign(at3=angle_df2.loc[:,'at3'].map(atom_id_remap))
    out_str = out_str + "\n\nAngles\n\n" + angle_df2[['type', 'at1', 'at2', 'at3']].to_string(header=False)
    
    dihedral_set1 = set(dihedral_df[dihedral_df['at1'].isin(selected_atoms)].index.to_list())
    dihedral_set2 = set(dihedral_df[dihedral_df['at2'].isin(selected_atoms)].index.to_list())
    dihedral_set3 = set(dihedral_df[dihedral_df['at3'].isin(selected_atoms)].index.to_list())
    dihedral_set4 = set(dihedral_df[dihedral_df['at4'].isin(selected_atoms)].index.to_list())
    dihedral_list = list(dihedral_set1 & dihedral_set2 & dihedral_set3 & dihedral_set4)
    dihedral_df2 = dihedral_df.loc[dihedral_list, :].reset_index(drop=True).copy(deep=True)
    dihedral_df2.index = dihedral_df2.index + 1
    dihedral_df2 = dihedral_df2.assign(at1=dihedral_df2.loc[:,'at1'].map(atom_id_remap))
    dihedral_df2 = dihedral_df2.assign(at2=dihedral_df2.loc[:,'at2'].map(atom_id_remap))
    dihedral_df2 = dihedral_df2.assign(at3=dihedral_df2.loc[:,'at3'].map(atom_id_remap))
    dihedral_df2 = dihedral_df2.assign(at4=dihedral_df2.loc[:,'at4'].map(atom_id_remap))
    out_str = out_str + "\n\nDihedrals\n\n" + dihedral_df2[['type', 'at1', 'at2', 'at3', 'at4']].to_string(header=False)
    
    improper_set1 = set(improper_df[improper_df['at1'].isin(selected_atoms)].index.to_list())
    improper_set2 = set(improper_df[improper_df['at2'].isin(selected_atoms)].index.to_list())
    improper_set3 = set(improper_df[improper_df['at3'].isin(selected_atoms)].index.to_list())
    improper_set4 = set(improper_df[improper_df['at4'].isin(selected_atoms)].index.to_list())
    improper_list = list(improper_set1 & improper_set2 & improper_set3 & improper_set4)
    improper_df2 = improper_df.loc[improper_list, :].reset_index(drop=True).copy(deep=True)
    improper_df2.index = improper_df2.index + 1
    improper_df2 = improper_df2.assign(at1=improper_df2.loc[:,'at1'].map(atom_id_remap))
    improper_df2 = improper_df2.assign(at2=improper_df2.loc[:,'at2'].map(atom_id_remap))
    improper_df2 = improper_df2.assign(at3=improper_df2.loc[:,'at3'].map(atom_id_remap))
    improper_df2 = improper_df2.assign(at4=improper_df2.loc[:,'at4'].map(atom_id_remap))
    out_str = out_str + "\n\nImpropers\n\n" + improper_df2[['type', 'at1', 'at2', 'at3', 'at4']].to_string(header=False)

    header_str = 'EC templated, generated by xyan11@uic.edu\n\n' + str(len(out_type)) + ' atoms\n' + \
    str(len(bond_df2)) + ' bonds\n' + \
    str(len(angle_df2)) + ' angles\n' + \
    str(len(dihedral_df2)) + ' dihedrals\n' + \
    str(len(improper_df2)) + ' impropers'

    sec_strs = [header_str,
               out_type.to_string(header=False), 
               out_q.to_string(header=False), 
               out_xyz.to_string(header=False), 
               bond_df2[['type', 'at1', 'at2']].to_string(header=False), 
               angle_df2[['type', 'at1', 'at2', 'at3']].to_string(header=False), 
               dihedral_df2[['type', 'at1', 'at2', 'at3', 'at4']].to_string(header=False), 
               improper_df2[['type', 'at1', 'at2', 'at3', 'at4']].to_string(header=False),
              ]
    sec_names = ['', 'Types', 'Charges', 'Coords', 'Bonds', 'Angles', 'Dihedrals', 'Impropers',]

    write_lmp_template(template_file_name, sec_names, sec_strs)
    return atom_id_remap

def build_map_file(map_fname, 
                   template_df, 
                   initiator_atom_ids, 
                   edge_atom_ids, 
                   delete_atom_template_ids, 
                   edge_atom_template_ids):
    mol1 = template_df[template_df["mol"]==1]
    mol2 = template_df[template_df["mol"]==2]
    initiator_atom_template_ids = template_df[template_df["id"].isin(initiator_atom_ids)]["template_atom_id"].tolist()
    outstr = """LAMMPS reaction map file generated by xyan11@uic.edu

""" + "%d" % len(edge_atom_template_ids) + """ edgeIDs
""" + "%d" % len(delete_atom_template_ids) + """ deleteIDs
""" + "%d" % len(template_df) + """ equivalences

InitiatorIDs

""" + "\n".join(["%d" % x for x in initiator_atom_template_ids]) + """

DeleteIDs

""" + "\n".join(["%d" % x for x in delete_atom_template_ids]) + """

EdgeIDs

""" + "\n".join(["%d" % x for x in edge_atom_template_ids]) + """

Equivalences

""" + "\n".join(["%d" % x + " " + "%d" % x for x in template_df["template_atom_id"]]) + """
"""
    with open(map_fname, "w") as wf:
        wf.write(outstr)

def build_peo(dop):
    monomer_mass_df, monomer_df, monomer_bond_df, monomer_angle_df, monomer_dihedral_df, monomer_improper_df, \
    [[monomer_xlo, monomer_xhi], [monomer_ylo, monomer_yhi], [monomer_zlo, monomer_zhi]], \
    monomer_lj_coeff, monomer_bond_coeff, monomer_angle_coeff, monomer_dihedral_coeff, monomer_improper_coeff = open_lmp_data('./eo_monomer.lmp', viz=False)
    monomer_cgenff = ['CG321', 'HGA2', 'HGA2', 'OG301', 'CG321', 'HGA2', 'HGA2', 'X']
    monomer_df['cgenff'] = monomer_cgenff

    head_mass_df, head_df, head_bond_df, head_angle_df, head_dihedral_df, head_improper_df, \
    [[head_xlo, head_xhi], [head_ylo, head_yhi], [head_zlo, head_zhi]], \
    head_lj_coeff, head_bond_coeff, head_angle_coeff, head_dihedral_coeff, head_improper_coeff = open_lmp_data('./MeOH_head.lmp', viz=False)
    head_cgenff = ['OG311', 'HGP1', 'CG321', 'HGA2', 'HGA2', 'X']
    head_df['cgenff'] = head_cgenff

    tail_mass_df, tail_df, tail_bond_df, tail_angle_df, tail_dihedral_df, tail_improper_df, \
    [[tail_xlo, tail_xhi], [tail_ylo, tail_yhi], [tail_zlo, tail_zhi]], \
    tail_lj_coeff, tail_bond_coeff, tail_angle_coeff, tail_dihedral_coeff, tail_improper_coeff = open_lmp_data('./MeOH_tail.lmp', viz=False)
    tail_cgenff = ['CG321', 'HGA2', 'HGA2', 'OG311', 'HGP1']
    tail_df['cgenff'] = tail_cgenff

    polymer_df = head_df.copy(deep=True)
    polymer_bond_df = head_bond_df.copy(deep=True)
    polymer_angle_df = head_angle_df.copy(deep=True)
    polymer_dihedral_df = head_dihedral_df.copy(deep=True)
    polymer_improper_df = head_improper_df.copy(deep=True)

    polymer_mass_df = head_mass_df.copy(deep=True)
    monomer_mass_df['type'] = monomer_mass_df['type'] + polymer_mass_df['type'].max() - 1
    monomer_df['type'] = monomer_df['type'] + polymer_mass_df['type'].max() - 1
    polymer_mass_df = pd.concat([polymer_mass_df.drop(polymer_mass_df.tail(1).index), monomer_mass_df], axis=0, ignore_index=True)
    polymer_mass_df.index = polymer_mass_df.index + 1

    tail_mass_df['type'] = tail_mass_df['type'] + polymer_mass_df['type'].max() - 1
    tail_df['type'] = tail_df['type'] + polymer_mass_df['type'].max() - 1
    polymer_mass_df = pd.concat([polymer_mass_df.drop(polymer_mass_df.tail(1).index), tail_mass_df], axis=0, ignore_index=True)
    polymer_mass_df.index = polymer_mass_df.index + 1

    for i in range(0, dop):
        tail_anchor = len(polymer_df)

        _df = monomer_df.copy(deep=True)
        head_atom = 1
        _bond_df = monomer_bond_df.copy(deep=True)
        _angle_df = monomer_angle_df.copy(deep=True)
        _dihedral_df = monomer_dihedral_df.copy(deep=True)
        _improper_df = monomer_improper_df.copy(deep=True)

        # process coordinates
        # reflection about X-Y plane
        if i % 2 == 1:
            _df.loc[:, 'z'] = 0.0 - _df.loc[:, 'z']
        # translate
        disp_vec = polymer_df.loc[tail_anchor, 'x':'z'] - _df.loc[1, 'x':'z']
        _df.loc[:, 'x':'z'] = _df.loc[:, 'x':'z'] + disp_vec

        # drop the anchor row, it will be replaced by the new head atom from _df
        polymer_df = polymer_df.drop(polymer_df.tail(1).index)

        # update ids in the new monomer
        max_atom_id = polymer_df['id'].max()
        max_bond_id = polymer_bond_df['id'].max()
        max_angle_id = polymer_angle_df['id'].max()
        max_dihedral_id = polymer_dihedral_df['id'].max()
        max_improper_id = polymer_improper_df['id'].max()

        _df['id'] = _df['id'] + max_atom_id
        _bond_df['at1'] = _bond_df['at1'] + max_atom_id
        _bond_df['at2'] = _bond_df['at2'] + max_atom_id
        _angle_df['at1'] = _angle_df['at1'] + max_atom_id
        _angle_df['at2'] = _angle_df['at2'] + max_atom_id
        _angle_df['at3'] = _angle_df['at3'] + max_atom_id
        _dihedral_df['at1'] = _dihedral_df['at1'] + max_atom_id
        _dihedral_df['at2'] = _dihedral_df['at2'] + max_atom_id
        _dihedral_df['at3'] = _dihedral_df['at3'] + max_atom_id
        _dihedral_df['at4'] = _dihedral_df['at4'] + max_atom_id
        _improper_df['at1'] = _improper_df['at1'] + max_atom_id
        _improper_df['at2'] = _improper_df['at2'] + max_atom_id
        _improper_df['at3'] = _improper_df['at3'] + max_atom_id
        _improper_df['at4'] = _improper_df['at4'] + max_atom_id

        head_atom = 1 + max_atom_id
        _bond_df['id'] = _bond_df['id'] + max_bond_id
        _angle_df['id'] = _angle_df['id'] + max_angle_id
        _dihedral_df['id'] = _dihedral_df['id'] + max_dihedral_id
        _improper_df['id'] = _improper_df['id'] + max_improper_id

        # combine polymer and new monomer
        polymer_df = pd.concat([polymer_df, _df], axis=0, ignore_index=True)
        polymer_df.index = polymer_df.index + 1
        polymer_bond_df = pd.concat([polymer_bond_df, _bond_df], axis=0, ignore_index=True)
        polymer_bond_df.index = polymer_bond_df.index + 1
        polymer_angle_df = pd.concat([polymer_angle_df, _angle_df], axis=0, ignore_index=True)
        polymer_angle_df.index = polymer_angle_df.index + 1
        polymer_dihedral_df = pd.concat([polymer_dihedral_df, _dihedral_df], axis=0, ignore_index=True)
        polymer_dihedral_df.index = polymer_dihedral_df.index + 1
        polymer_improper_df = pd.concat([polymer_improper_df, _improper_df], axis=0, ignore_index=True)
        polymer_improper_df.index = polymer_improper_df.index + 1


    tail_anchor = len(polymer_df)

    _df = tail_df.copy(deep=True)
    head_atom = 1
    _bond_df = tail_bond_df.copy(deep=True)
    _angle_df = tail_angle_df.copy(deep=True)
    _dihedral_df = tail_dihedral_df.copy(deep=True)
    _improper_df = tail_improper_df.copy(deep=True)

    # process coordinates
    # reflection about X-Y plane
    if i % 2 != 1:
        _df.loc[:, 'z'] = 0.0 - _df.loc[:, 'z']
    # translate
    disp_vec = polymer_df.loc[tail_anchor, 'x':'z'] - _df.loc[1, 'x':'z']
    _df.loc[:, 'x':'z'] = _df.loc[:, 'x':'z'] + disp_vec

    # drop the anchor row, it will be replaced by the new head atom from _df
    polymer_df = polymer_df.drop(polymer_df.tail(1).index)

    # update ids in the new monomer
    max_atom_id = polymer_df['id'].max()
    max_bond_id = polymer_bond_df['id'].max()
    max_angle_id = polymer_angle_df['id'].max()
    max_dihedral_id = polymer_dihedral_df['id'].max()
    max_improper_id = polymer_improper_df['id'].max()

    _df['id'] = _df['id'] + max_atom_id
    _bond_df['at1'] = _bond_df['at1'] + max_atom_id
    _bond_df['at2'] = _bond_df['at2'] + max_atom_id
    _angle_df['at1'] = _angle_df['at1'] + max_atom_id
    _angle_df['at2'] = _angle_df['at2'] + max_atom_id
    _angle_df['at3'] = _angle_df['at3'] + max_atom_id
    _dihedral_df['at1'] = _dihedral_df['at1'] + max_atom_id
    _dihedral_df['at2'] = _dihedral_df['at2'] + max_atom_id
    _dihedral_df['at3'] = _dihedral_df['at3'] + max_atom_id
    _dihedral_df['at4'] = _dihedral_df['at4'] + max_atom_id
    _improper_df['at1'] = _improper_df['at1'] + max_atom_id
    _improper_df['at2'] = _improper_df['at2'] + max_atom_id
    _improper_df['at3'] = _improper_df['at3'] + max_atom_id
    _improper_df['at4'] = _improper_df['at4'] + max_atom_id

    head_atom = 1 + max_atom_id
    _bond_df['id'] = _bond_df['id'] + max_bond_id
    _angle_df['id'] = _angle_df['id'] + max_angle_id
    _dihedral_df['id'] = _dihedral_df['id'] + max_dihedral_id
    _improper_df['id'] = _improper_df['id'] + max_improper_id

    # combine polymer and new monomer
    polymer_df = pd.concat([polymer_df, _df], axis=0, ignore_index=True)
    polymer_df.index = polymer_df.index + 1
    polymer_bond_df = pd.concat([polymer_bond_df, _bond_df], axis=0, ignore_index=True)
    polymer_bond_df.index = polymer_bond_df.index + 1
    polymer_angle_df = pd.concat([polymer_angle_df, _angle_df], axis=0, ignore_index=True)
    polymer_angle_df.index = polymer_angle_df.index + 1
    polymer_dihedral_df = pd.concat([polymer_dihedral_df, _dihedral_df], axis=0, ignore_index=True)
    polymer_dihedral_df.index = polymer_dihedral_df.index + 1
    polymer_improper_df = pd.concat([polymer_improper_df, _improper_df], axis=0, ignore_index=True)
    polymer_improper_df.index = polymer_improper_df.index + 1

    return polymer_mass_df, polymer_df, polymer_bond_df, polymer_angle_df, polymer_dihedral_df, polymer_improper_df

def assign_cgenff(data_fname, polymer_mass_df, polymer_df, polymer_bond_df, polymer_angle_df, polymer_dihedral_df, polymer_improper_df):
    polymer_df['comment'] = '# ' + polymer_df['cgenff']
    id2cgenff = polymer_df.set_index('id').to_dict()['cgenff']
    polymer_bond_df['cgenff'] = polymer_bond_df['at1'].map(id2cgenff) + '-' + polymer_bond_df['at2'].map(id2cgenff)
    polymer_bond_df['comment'] = '# ' + polymer_bond_df['cgenff']
    polymer_angle_df['cgenff'] = polymer_angle_df['at1'].map(id2cgenff) + '-' + polymer_angle_df['at2'].map(id2cgenff) + '-' + polymer_angle_df['at3'].map(id2cgenff)
    polymer_angle_df['comment'] = '# ' + polymer_angle_df['cgenff']
    polymer_dihedral_df['cgenff'] = polymer_dihedral_df['at1'].map(id2cgenff) + '-' + polymer_dihedral_df['at2'].map(id2cgenff) + '-' + polymer_dihedral_df['at3'].map(id2cgenff) + '-' + polymer_dihedral_df['at4'].map(id2cgenff)
    polymer_dihedral_df['comment'] = '# ' + polymer_dihedral_df['cgenff']
    polymer_improper_df['cgenff'] = polymer_improper_df['at1'].map(id2cgenff) + '-' + polymer_improper_df['at2'].map(id2cgenff) + '-' + polymer_improper_df['at3'].map(id2cgenff) + '-' + polymer_improper_df['at4'].map(id2cgenff)
    polymer_improper_df['comment'] = '# ' + polymer_improper_df['cgenff']

    prm_str = None
    with open('./par_all36_cgenff.prm', 'r') as prm:
        prm_str = prm.read()

    lj_str = [x for x in prm_str.split('\nNONBONDED')[1].split('\nNBFIX')[0].split('!hydrogens')[1].split('\n') if (x and not x.startswith('!'))]
    lj_prm_list, lj_prm_comment = list(map(list, zip(*[x.split('!') for x in lj_str])))

    lj_prm = pd.read_csv(io.StringIO('\n'.join(lj_prm_list)),
                        comment='!', sep=r'\s+', usecols=[0, 1, 2, 3, 4, 5, 6],
                        names=['at', 'ignored1', 'epsilon', 'Rmin/2', 'ignored2', 'eps,1-4', 'Rmin/2,1-4']).fillna(0)
    # V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
    # epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
    # Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
    lj_prm['comment'] = lj_prm_comment
    lj_prm['comment'] = '# ' + lj_prm['comment']

    lj_prm = lj_prm[lj_prm['at'].isin(polymer_df['cgenff'])].reset_index(drop=True)
    lj_prm.index = lj_prm.index + 1
    lj_prm['sigma/2'] = lj_prm['Rmin/2'] / (2.**(1./6.))
    lj_prm['sigma/2,1-4'] = lj_prm['Rmin/2,1-4'] / (2.**(1./6.))


    bond_str = [x for x in prm_str.split('\nBONDS')[1].split('\nANGLES')[0].split('\n') if (x and not x.startswith('!'))]
    bond_prm_list, bond_prm_comment = list(map(list, zip(*[x.split('!') for x in bond_str])))

    bond_prm = pd.read_csv(io.StringIO('\n'.join(bond_prm_list)),
                        comment='!', sep=r'\s+', usecols=[0, 1, 2, 3],
                        names=['at1', 'at2', 'Kb', 'b0'])
    # V(bond) = Kb(b - b0)**2
    # Kb: kcal/mole/A**2
    # b0: A
    bond_prm['comment'] = bond_prm_comment
    bond_prm['comment'] = '# ' + bond_prm['comment']

    bond_prm1 = bond_prm[bond_prm['at1'].isin(polymer_df['cgenff'])].index
    bond_prm2 = bond_prm[bond_prm['at2'].isin(polymer_df['cgenff'])].index

    bond_index = bond_prm1.intersection(bond_prm2)
    bond_prm = bond_prm.loc[bond_index].reset_index(drop=True)
    bond_prm.index = bond_prm.index + 1

    angle_str = [x for x in prm_str.split('\nANGLES')[1].split('\nDIHEDRALS')[0].split('\n') if (x and not x.startswith('!'))]
    angle_prm_list, angle_prm_comment = list(map(list, zip(*[x.split('!') for x in angle_str])))

    angle_prm = pd.read_csv(io.StringIO('\n'.join(angle_prm_list)),
                            comment='!', sep=r'\s+', usecols=[0, 1, 2, 3, 4, 5, 6],
                            names=['at1', 'at2', 'at3', 'Ktheta', 'r', 'kUB', 'rUB']).fillna(0)
    #V(angle) = Ktheta(Theta - Theta0)**2
    #Ktheta: kcal/mole/rad**2
    #Theta0: degrees
    #V(Urey-Bradley) = Kub(S - S0)**2
    #Kub: kcal/mole/A**2 (Urey-Bradley)
    angle_prm['comment'] = angle_prm_comment
    angle_prm['comment'] = '# ' + angle_prm['comment']

    angle_prm1 = angle_prm[angle_prm['at1'].isin(polymer_df['cgenff'])].index
    angle_prm2 = angle_prm[angle_prm['at2'].isin(polymer_df['cgenff'])].index
    angle_prm3 = angle_prm[angle_prm['at3'].isin(polymer_df['cgenff'])].index
    ang_index = angle_prm1.intersection(angle_prm2).intersection(angle_prm3)
    angle_prm = angle_prm.loc[ang_index].reset_index(drop=True)
    angle_prm.index = angle_prm.index + 1

    dihedral_str = [x for x in prm_str.split('\nDIHEDRALS')[1].split('\nIMPROPERS')[0].split('\n') if (x and not x.startswith('!'))]
    dihedral_prm_list, dihedral_prm_comment = list(map(list, zip(*[x.split('!') for x in dihedral_str])))
    dihedral_prm = pd.read_csv(io.StringIO('\n'.join(dihedral_prm_list)),
                            comment='!', sep=r'\s+', usecols=[0, 1, 2, 3, 4, 5, 6],
                            names=['at1', 'at2', 'at3', 'at4', 'Kchi', 'n', 'delta'])
    # V(dihedral) = Kchi(1 + cos(n(chi) - delta))
    # Kchi: kcal/mole
    # n: multiplicity
    # delta: degrees
    dihedral_prm['comment'] = dihedral_prm_comment
    dihedral_prm['comment'] = '# ' + dihedral_prm['comment']

    dihedral_prm1 = dihedral_prm[dihedral_prm['at1'].isin(polymer_df['cgenff'])].index
    dihedral_prm2 = dihedral_prm[dihedral_prm['at2'].isin(polymer_df['cgenff'])].index
    dihedral_prm3 = dihedral_prm[dihedral_prm['at3'].isin(polymer_df['cgenff'])].index
    dihedral_prm4 = dihedral_prm[dihedral_prm['at4'].isin(polymer_df['cgenff'])].index
    dih_index = dihedral_prm1.intersection(dihedral_prm2).intersection(dihedral_prm3).intersection(dihedral_prm4)
    dihedral_prm = dihedral_prm.loc[dih_index].reset_index(drop=True)
    dihedral_prm.index = dihedral_prm.index + 1

    improper_str = [x for x in prm_str.split('\nIMPROPERS')[1].split('\nNONBONDED')[0].split('\n') if (x and not x.startswith('!'))]
    improper_prm_list, improper_prm_comment = list(map(list, zip(*[x.split('!') for x in improper_str])))
    improper_prm = pd.read_csv(io.StringIO('\n'.join(improper_prm_list)),
                            comment='!', sep=r'\s+', usecols=[0, 1, 2, 3, 4, 5, 6],
                            names=['at1', 'at2', 'at3', 'at4', 'Kpsi', 'ignored', 'psi0'])
    # V(improper) = Kpsi(psi - psi0)**2
    # Kpsi: kcal/mole/rad**2
    # psi0: degrees
    # note that the second column of numbers (0) is ignored
    improper_prm['comment'] = improper_prm_comment
    improper_prm['comment'] = '# ' + improper_prm['comment']

    improper_prm1 = improper_prm[improper_prm['at1'].isin(polymer_df['cgenff'])].index
    improper_prm2 = improper_prm[improper_prm['at2'].isin(polymer_df['cgenff'])].index
    improper_prm3 = improper_prm[improper_prm['at3'].isin(polymer_df['cgenff'])].index
    improper_prm4 = improper_prm[improper_prm['at4'].isin(polymer_df['cgenff'])].index
    imp_index = improper_prm1.intersection(improper_prm2).intersection(improper_prm3).intersection(improper_prm4)
    improper_prm = improper_prm.loc[imp_index].reset_index(drop=True)
    improper_prm.index = improper_prm.index + 1

    self_lj = []
    _polymer_df = polymer_df.drop_duplicates(subset=['type']).reset_index(drop=True)
    _polymer_df.index = _polymer_df.index + 1
    for i in range(1, 1+len(_polymer_df)):
        atom_type = _polymer_df.loc[i, 'cgenff']
        curr_lj_prm = lj_prm[lj_prm['at']==atom_type]
        self_lj.append(curr_lj_prm)
    self_lj_df = pd.concat(self_lj, axis=0).reset_index(drop=True)
    self_lj_df.index = self_lj_df.index + 1
    self_lj_df['id'] = self_lj_df.index

    pair_coeff = ''
    for i in range(1, 1+len(_polymer_df)):
        for j in range(i, 1+len(_polymer_df)):
            mixed_sigma = self_lj_df.at[i, 'sigma/2'] + self_lj_df.at[j, 'sigma/2']
            mixed_sigma14 = self_lj_df.at[i, 'sigma/2,1-4'] + self_lj_df.at[j, 'sigma/2,1-4']
            mixed_epsilon = (self_lj_df.at[i, 'epsilon'] * self_lj_df.at[j, 'epsilon']) ** (1./2.)
            mixed_epsilon14 = (self_lj_df.at[i, 'eps,1-4'] * self_lj_df.at[j, 'eps,1-4']) ** (1./2.)

            pair_coeff = pair_coeff + 'pair_coeff ' + "{:.0f}".format(self_lj_df.at[i, 'id']) + ' ' + \
                                                    "{:.0f}".format(self_lj_df.at[j, 'id']) + ' ' + \
                                                    "{:.6f}".format(mixed_epsilon) + ' ' + \
                                                    "{:.6f}".format(mixed_sigma) + ' ' + \
                                                    "{:.6f}".format(mixed_epsilon) + ' ' + \
                                                    "{:.6f}".format(mixed_epsilon14) + ' ' + \
                                                    '# ' + self_lj_df.at[i, 'at'] + '-' + self_lj_df.at[j, 'at'] + '\n'
    #print(pair_coeff)
    _polymer_bond_df = polymer_bond_df.drop_duplicates(subset=['type']).reset_index(drop=True)
    _polymer_bond_df.index = _polymer_bond_df.index + 1
    bond_coeff = ""
    for i in range(1, 1+len(_polymer_bond_df)):
        all_atoms = _polymer_bond_df.loc[i, 'cgenff'].split('-')
        prm_found = False
        for j in range(1, 1+len(bond_prm)):
            if bond_prm.loc[j, ['at1', 'at2']].to_list() == all_atoms or bond_prm.loc[j, ['at2', 'at1']].to_list() == all_atoms:
                bond_coeff = bond_coeff + 'bond_coeff ' + str(i) + ' ' + \
                                "{:.4f}".format(bond_prm.at[j, 'Kb']) + ' ' + \
                                "{:.4f}".format(bond_prm.at[j, 'b0']) + ' ' + \
                                bond_prm.at[j, 'comment'] + '\n'
                prm_found = True
                break
        if prm_found == False:
            bond_coeff = bond_coeff + 'bond_coeff ' + str(i) + ' none' + '\n'
    #print(bond_coeff)

    angle_coeff = ""
    _polymer_angle_df = polymer_angle_df.drop_duplicates(subset=['type']).reset_index(drop=True)
    _polymer_angle_df.index = _polymer_angle_df.index + 1
    for i in range(1, 1+len(_polymer_angle_df)):
        all_atoms = _polymer_angle_df.loc[i, 'cgenff'].split('-')
        prm_found = False
        for j in range(1, 1+len(angle_prm)):
            if angle_prm.loc[j, ['at1', 'at2', 'at3']].to_list() == all_atoms or angle_prm.loc[j, ['at3', 'at2', 'at1']].to_list() == all_atoms:
                angle_coeff = angle_coeff + 'angle_coeff ' + str(i) + ' ' + \
                                "{:.4f}".format(angle_prm.at[j, 'Ktheta']) + ' ' + \
                                "{:.2f}".format(angle_prm.at[j, 'r']) + ' ' + \
                                "{:.4f}".format(angle_prm.at[j, 'kUB']) + ' ' + \
                                "{:.4f}".format(angle_prm.at[j, 'rUB']) + ' ' + \
                                angle_prm.at[j, 'comment'] + '\n'
                prm_found = True
                break
        if prm_found == False:
            angle_coeff = angle_coeff + 'angle_coeff ' + str(i) + ' none' + '\n'
    #print(angle_coeff)

    dihedral_coeff = ""
    _polymer_dihedral_df = polymer_dihedral_df.drop_duplicates(subset=['type']).reset_index(drop=True)
    _polymer_dihedral_df.index = _polymer_dihedral_df.index + 1
    for i in range(1, 1+len(_polymer_dihedral_df)):
        all_atoms = _polymer_dihedral_df.loc[i, 'cgenff'].split('-')
        prm_found = False
        for j in range(1, 1+len(dihedral_prm)):
            if dihedral_prm.loc[j, ['at1', 'at2', 'at3', 'at4']].to_list() == all_atoms or dihedral_prm.loc[j, ['at4', 'at3', 'at2', 'at1']].to_list() == all_atoms:
                weight = 1.0
                if all_atoms[0].startswith('C') and all_atoms[1].startswith('C') and all_atoms[2].startswith('C') and all_atoms[3].startswith('C'):
                    weight = 0.5
                dihedral_coeff = dihedral_coeff + 'dihedral_coeff ' + str(i) + ' ' + \
                                "{:.4f}".format(dihedral_prm.at[j, 'Kchi']) + ' ' + \
                                "{:.0f}".format(dihedral_prm.at[j, 'n']) + ' ' + \
                                "{:.0f}".format(dihedral_prm.at[j, 'delta']) + ' ' + \
                                str(weight) + ' ' + \
                                dihedral_prm.at[j, 'comment'] + '\n'
                prm_found = True
                break
        if prm_found == False:
            dihedral_coeff = dihedral_coeff + 'dihedral_coeff ' + str(i) + ' none' + '\n'
    #print(dihedral_coeff)

    improper_coeff = ""
    _polymer_improper_df = polymer_improper_df.drop_duplicates(subset=['type']).reset_index(drop=True)
    _polymer_improper_df.index = _polymer_improper_df.index + 1
    for i in range(1, 1+len(_polymer_improper_df)):
        all_atoms = _polymer_improper_df.loc[i, 'cgenff'].split('-')
        center_atom = all_atoms[0]
        satellite_atoms = sorted(all_atoms[1:4])
        search_df = improper_prm[improper_prm['at1'] == center_atom].reset_index(drop=True)
        if len(search_df) > 0:
            for j in range(0, len(search_df)):
                satellite_atoms_prm = sorted(search_df.loc[j, ['at2', 'at3', 'at4']])
                if satellite_atoms == satellite_atoms_prm:
                    #if [all_atoms[1], all_atoms[2], all_atoms[3]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        #print('No swap!')
                    if [all_atoms[1], all_atoms[3], all_atoms[2]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        improper_df.loc[i, ['at2','at3','at4']] = improper_df.loc[i, ['at2','at4','at3']]
                    elif [all_atoms[2], all_atoms[1], all_atoms[3]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        improper_df.loc[i, ['at2','at3','at4']] = improper_df.loc[i, ['at3','at2','at4']]
                    elif [all_atoms[3], all_atoms[2], all_atoms[1]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        improper_df.loc[i, ['at2','at3','at4']] = improper_df.loc[i, ['at4','at3','at2']]
                    elif [all_atoms[2], all_atoms[3], all_atoms[1]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        improper_df.loc[i, ['at2','at3','at4']] = improper_df.loc[i, ['at3','at4','at2']]
                    elif [all_atoms[3], all_atoms[1], all_atoms[2]] == search_df.loc[j, ['at2', 'at3', 'at4']].to_list():
                        improper_df.loc[i, ['at2','at3','at4']] = improper_df.loc[i, ['at4','at2','at3']]
                    improper_df.loc[i, 'cgenff'] = id2cgenff[improper_df.loc[i, 'at1']] + '-' + \
                                                id2cgenff[improper_df.loc[i, 'at2']] + '-' + \
                                                id2cgenff[improper_df.loc[i, 'at3']] + '-' + \
                                                id2cgenff[improper_df.loc[i, 'at4']]
                    improper_df.loc[i,'comment'] = '# ' + improper_df.loc[i,'cgenff']
                    improper_coeff = improper_coeff + 'improper_coeff ' + str(i) + ' ' + \
                                    "{:.4f}".format(search_df.at[j, 'Kpsi']) + ' ' + \
                                    "{:.2f}".format(search_df.at[j, 'psi0']) + ' ' + \
                                    search_df.at[j, 'comment'] + '\n'
                    break
        else:
            improper_coeff = improper_coeff + 'improper_coeff ' + str(i) + ' 0.0000 0.00' + '\n'


    cgenff36 = '''pair_style     lj/charmmfsw/coul/long 10.0 12.0
pair_modify    tail yes
kspace_style   pppm 1.0e-5
bond_style     harmonic
angle_style    charmm
dihedral_style charmmfsw
improper_style harmonic
special_bonds  charmm
neighbor       2.0 bin
neigh_modify   every 1 delay 0 check yes
    '''

    cgenff36 = cgenff36 + '''
    ''' + pair_coeff + '''

    ''' + bond_coeff + '''

    ''' + angle_coeff + '''

    ''' + dihedral_coeff + '''

    ''' + improper_coeff + '''

    '''

    with open('./cgenff36.ff', 'w') as ff_file:
        ff_file.write(cgenff36)

    min_coord = min([polymer_df['x'].min() - 1.0, polymer_df['y'].min() - 1.0, polymer_df['z'].min() - 1.0])
    max_coord = max([polymer_df['x'].max() + 1.0, polymer_df['y'].max() + 1.0, polymer_df['z'].max() + 1.0])

    xlo = polymer_df['x'].min() - 1.0
    xhi = polymer_df['x'].max() + 1.0
    ylo = polymer_df['y'].min() - 1.0
    yhi = polymer_df['y'].max() + 1.0
    zlo = polymer_df['z'].min() - 1.0
    zhi = polymer_df['z'].max() + 1.0

    new_sec_strs = ['LAMMPS data file for cgenff36 by xyan11@uic.edu' + """

    """ + str(len(polymer_df)) + """ atoms
    """ + str(len(polymer_mass_df)) + """ atom types
    """ + str(len(polymer_bond_df)) + """ bonds
    """ + str(int(polymer_bond_df['type'].max())) + """ bond types
    """ + str(len(polymer_angle_df)) + """ angles
    """ + str(int(polymer_angle_df['type'].max())) + """ angle types
    """ + str(len(polymer_dihedral_df)) + """ dihedrals
    """ + str(int(polymer_dihedral_df['type'].max())) + """ dihedral types
    """ + str(len(polymer_improper_df)) + """ impropers
    """ + str(int(polymer_improper_df['type'].max())) + """ improper types

    """ + "{:.6f}".format(xlo) + """ """ + "{:.6f}".format(xhi) + """ xlo xhi
    """ + "{:.6f}".format(ylo) + """ """ + "{:.6f}".format(yhi) + """ ylo yhi
    """ + "{:.6f}".format(zlo) + """ """ + "{:.6f}".format(zhi) + """ zlo zhi
    """,
    polymer_mass_df[['type', 'mass']].to_string(header=None, index=None) + '\n',
    polymer_df[['id', 'mol', 'type', 'q', 'x', 'y', 'z', 'comment']].to_string(header=None, index=None) + '\n',
    polymer_bond_df[['id', 'type', 'at1', 'at2', 'comment']].to_string(header=None, index=None) + '\n',
    polymer_angle_df[['id', 'type', 'at1', 'at2', 'at3', 'comment']].to_string(header=None, index=None) + '\n',
    polymer_dihedral_df[['id', 'type', 'at1', 'at2', 'at3', 'at4', 'comment']].to_string(header=None, index=None) + '\n',
    polymer_improper_df[['id', 'type', 'at1', 'at2', 'at3', 'at4', 'comment']].to_string(header=None, index=None) + '\n']

    new_sec_names = ['head',
    '\n\nMasses\n\n',
    '\n\nAtoms \n\n',
    '\n\nBonds \n\n',
    '\n\nAngles \n\n',
    '\n\nDihedrals\n\n',
    '\n\nImpropers\n\n']
    write_lmp_data(data_fname, new_sec_names, new_sec_strs)

    return lj_prm, bond_prm, angle_prm, dihedral_prm, improper_prm, xlo, xhi, ylo, yhi, zlo, zhi



def map_discrete_colors(ncolors):
    n_per_channel = int((ncolors - 1) / 3.) + 1
    R = (np.concatenate([np.zeros(n_per_channel), 
                         np.linspace(0, 1, n_per_channel),
                         np.ones(n_per_channel)], axis=0) * 255).astype(int)
    G = (np.concatenate([np.linspace(0, 1, n_per_channel), 
                         np.ones(n_per_channel), 
                         np.linspace(1, 0, n_per_channel)], axis=0) * 255).astype(int)
    B = (np.concatenate([np.ones(n_per_channel), 
                         np.linspace(1, 0, n_per_channel),
                         np.zeros(n_per_channel)], axis=0) * 255).astype(int)
    A = np.ones(n_per_channel * 3).astype(int)
    colors = ["rgba(" + ",".join(x) + ")" for x in (np.array([R, G, B, A]).T).astype(str).tolist()]
    return dict(zip([x+1 for x in range(ncolors)], colors))

def read_lmp_fix_file(fname):
    outlines = []
    with io.open(fname, "r", newline="\n") as rf:
        outlines = rf.readlines()
    line2_args = [x.strip() for x in list(filter(None, outlines[1].replace("#", "").split(" ")))]
    
    args = dict(zip(line2_args, [[] for x in line2_args]))
    time_arg_name = line2_args[0]
    nlines_arg_name = line2_args[1]
    conserve_arg_name = line2_args[2]
    args["dfs"] = []
    
    df_header = [x.strip() for x in list(filter(None, outlines[2].replace("#", "").split(" ")))]

    line_i = 3
    while line_i < len(outlines):
        new_args = [int(np.round(float(x), 0)) for x in list(filter(None, outlines[line_i].replace("#", "").split(" ")))]
        args[time_arg_name].append(new_args[0])
        args[nlines_arg_name].append(new_args[1])
        args[conserve_arg_name].append(new_args[2])
        df = pd.read_csv(io.StringIO("".join(outlines[line_i+1:line_i+1+new_args[1]])), sep=r"\s+", names=df_header)
        args["dfs"].append(df)
        line_i = line_i + new_args[1] + 1
    return args