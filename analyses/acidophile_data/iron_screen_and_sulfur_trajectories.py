"""Iron-screen classifications and soluble-sulfur endpoint calculations."""
from __future__ import annotations
import json
import re
import pandas as pd
from .plate_comparisons import CONDITION_KEYS
from .peak_response import plate_means
from .publication_response_validation import consecutive_snr_summary


SULFATE_CATIONS = {
    'aluminum':('Al',3),'cadmium':('Cd',2),'cerium_iii':('Ce',3),'cerium_iv':('Ce',4),
    'chromium_iii':('Cr',3),'cobalt_ii':('Co',2),'copper_ii':('Cu',2),
    'dysprosium_iii':('Dy',3),'erbium_iii':('Er',3),'europium_iii':('Eu',3),
    'gadolinium_iii':('Gd',3),'gallium':('Ga',3),'indium':('In',3),'iron_ii':('Fe',2),
    'iron_iii':('Fe',3),'lanthanum_iii':('La',3),'lithium':('Li',1),'magnesium':('Mg',2),
    'manganese_ii':('Mn',2),'neodymium_iii':('Nd',3),'nickel_ii':('Ni',2),
    'praseodymium_iii':('Pr',3),'samarium_iii':('Sm',3),'silver':('Ag',1),
    'terbium_iii':('Tb',3),'ytterbium_iii':('Yb',3),'yttrium_iii':('Y',3),
    'zinc':('Zn',2),'zirconium_iv':('Zr',4),'sodium':('Na',1),'potassium':('K',1),
}


FORMULAS = {
    'ammonium_sulfate':'(NH4)2SO4','iron_iii_chloride':'FeCl3',
    'monopotassium_phosphate':'KH2PO4','potassium_bromide':'KBr','potassium_chloride':'KCl',
    'potassium_dichromate':'K2Cr2O7','potassium_fluoride':'KF','sodium_gluconate':'NaC6H11O7',
    'potassium_iodate':'KIO3','potassium_iodide':'KI','potassium_nitrate':'KNO3',
    'potassium_perrhenate':'KReO4','potassium_tetrathionate':'K2S4O6',
    'potassium_thiocyanate':'KSCN','sodium_chloride':'NaCl','sodium_chromate':'Na2CrO4',
    'sodium_formate':'HCOONa','sodium_molybdate':'Na2MoO4','sodium_selenite':'Na2SeO3',
    'sodium_sulfite':'Na2SO3','sodium_tungstate':'Na2WO4','sucrose':'C12H22O11',
    'titanyl_sulfate':'TiOSO4','vanadyl_sulfate':'VOSO4','citric_acid':'C6H8O7',
    'gallium':'Ga','sodium_thiosulfate':'Na2S2O3','sulfuric_acid':'H2SO4',
}


for key,(symbol,charge) in SULFATE_CATIONS.items():
    FORMULAS[key+'_sulfate'] = (symbol+'2SO4' if charge==1 else symbol+'SO4' if charge==2
                              else symbol+'2(SO4)3' if charge==3 else symbol+'(SO4)2')


def chemical(component):
    formula=FORMULAS.get(component)
    if formula is None: return component.replace('_',' ')
    return r'$\mathrm{'+re.sub(r'(\d+)',r'_{\1}',formula)+'}$'


def row_label(row):
    component=row.get('principal_component','')
    quantities=json.loads(row.get('component_final_mM_json','{}'))
    if not component:component='citric_acid' if 'citric_acid' in quantities else next(iter(quantities),'control')
    label=chemical(component)
    if 'citric_acid' in quantities and 'gallium' in quantities:
        label += f' / {quantities["gallium"]:g} mM Ga'
    ph=row.get('nominal_initial_pH_json','null')
    if pd.notna(ph) and str(ph) not in ('null','nan',''):
        label+=f' / pH {ph}'
    return label


def screen(tables):
    rows=[]
    frames=[tables['af_iron_plate_trajectories'],tables['validation_lf_plate_trajectories']]
    for f in frames:
        for identity,g in f[f.signal.eq('raw_A380')].groupby(CONDITION_KEYS+['batch_id','plate_id']):
            g=g.sort_values('days_mean');finite=g[g.biotic_change.notna()]
            maxrow=finite.loc[finite.biotic_change.idxmax()] if len(finite) else g.iloc[0]
            event=g.biotic_change.ge(1)&g.abiotic_change.lt(1)
            r=dict(zip(CONDITION_KEYS+['batch_id','plate_id'],identity)) | dict(
                maximum_biotic_change=maxrow.biotic_change,abiotic_at_maximum=maxrow.abiotic_change,
                net_at_maximum=maxrow.net_change,maximum_time_d=maxrow.days_mean,
                optical_activity_event=bool(event.any()),candidate_acquisitions=len(g),
                eligible_biotic_acquisitions=g.biotic_change.count(),
                eligible_matched_acquisitions=int((g.biotic_change.notna()&g.abiotic_change.notna()).sum()),
                culture_ids=maxrow.culture_ids)
            for field in ['component_final_mM_json','medium_id','condition_json']:
                if field in g:r[field]=g[field].iloc[0]
            rows.append(r)
    f=pd.DataFrame(rows)
    af=tables['af_iron_condition_coverage']
    metadata=af[CONDITION_KEYS+['principal_component','principal_final_mM','additive_class',
        'component_final_mM_json','nominal_initial_pH_json','medium_id']].copy()
    lf=tables['validation_lf_plate_trajectories'][CONDITION_KEYS+['component_final_mM_json','medium_id']].drop_duplicates()
    lf['principal_component']=lf.component_final_mM_json.map(lambda s:next(iter(json.loads(s))))
    lf['principal_final_mM']=lf.apply(lambda r:json.loads(r.component_final_mM_json)[r.principal_component],axis=1)
    lf['additive_class']='LF additions';lf['nominal_initial_pH_json']='null'
    metadata=pd.concat([metadata,lf],ignore_index=True)
    f=f.drop(columns=['component_final_mM_json','medium_id','condition_json'],errors='ignore').merge(metadata,on=CONDITION_KEYS,validate='many_to_one')
    f['chemical_label']=f.apply(row_label,axis=1);f['dose_M']=f.principal_final_mM/1000
    group=CONDITION_KEYS+['chemical_label','dose_M','additive_class','medium_id']
    summary=f.groupby(group,as_index=False).agg(mean=('maximum_biotic_change','mean'),
        sd=('maximum_biotic_change','std'),eligible_plate_n=('maximum_biotic_change','count'),
        candidate_plate_n=('plate_id','size'),activity_event_n=('optical_activity_event','sum'),
        batch_n=('batch_id','nunique'))
    return {'screen_plate_evidence':f,'screen_summary':summary}


def sulfur(tables):
    source=tables['validation_sulfur_trajectory_evidence'];curves=plate_means(source)
    keys=CONDITION_KEYS+['batch_id','plate_id','signal','ordinal']
    definitions=['signed','magnitude','signed_common','magnitude_common']
    spec={f'technical_snr_{name}':(f'technical_snr_{name}','mean') for name in definitions}
    spec.update(technical_snr_n=('technical_snr_signed','count'),
                technical_snr_common_n=('technical_snr_signed_common','count'))
    signal=source.groupby(keys,as_index=False).agg(**spec)
    signal['technical_snr']=signal.technical_snr_signed
    curves=curves.merge(signal,on=keys,validate='one_to_one')
    metadata=tables['context_sulfur_eligibility'][CONDITION_KEYS+['component_final_mM_json','medium_id']].drop_duplicates()
    curves=curves.merge(metadata,on=CONDITION_KEYS,validate='many_to_one')
    rows=[]
    for identity,g in curves.groupby(CONDITION_KEYS+['batch_id','plate_id','signal']):
        g=g.sort_values('ordinal');last=g.iloc[-1]
        q=json.loads(last.component_final_mM_json)
        component=next((c for c in q if '_'+c+'_OD' in last.experiment_id),next(iter(q)))
        result=dict(zip(CONDITION_KEYS+['batch_id','plate_id','signal'],identity))
        for name in definitions:
            result[f'sustained_snr_{name}'],result[f'consecutive_pair_n_{name}']=(
                consecutive_snr_summary(g,f'technical_snr_{name}'))
            result[f'final_snr_{name}']=last[f'technical_snr_{name}']
        result.update(
            dict(sustained_snr=result['sustained_snr_signed'],
                final_biotic_change=last.biotic_change,final_abiotic_change=last.abiotic_change,
                final_net_change=last.net_change,final_day=last.days_mean,
                final_ordinal=last.ordinal,final_snr=result['final_snr_signed'],
                final_snr_n=int(last.technical_snr_n),
                final_snr_common_n=int(last.technical_snr_common_n),
                expected_well_n=int(last.expected_well_n),culture_ids=last.culture_ids,
                chemical_label=chemical(component),dose_M=q[component]/1000,
                medium_id=last.medium_id,component_final_mM_json=last.component_final_mM_json))
        rows.append(result)
    return {'sulfur_plate_evidence':pd.DataFrame(rows),
            'sulfur_plate_trajectories':curves}
