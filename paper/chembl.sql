select h1,count(*) as cnt from chembl28.spectral_hk group by h1 order by cnt DESC
;

select h1,count(*) as cnt from chembl28_ncats.spectral_hk group by h1 order by cnt DESC
;

select d.assay_id,d.description,type,standard_units,count(distinct b.molregno) as N, avg(log(10,b.standard_value*10^-9)), std(log(10,b.standard_value))
from chembl_id_lookup a, activities b, chembl28_ncats.spectral_hk c,assays d 
where a.entity_id = b.molregno
and a.chembl_id = c.chembl_id
and c.h1='KKPV1D97H'
#and c.h2='NN51MHRR4DXYNAQ6H78'
#and c.h2='NKVWAXFWSFGCG8MFJLT'
and standard_relation='='
and b.assay_id=d.assay_id
group by d.assay_id,d.description,type,units
having count(*) > 1
order by N desc
;

select d.description,e.canonical_smiles,type,standard_units,b.standard_value
from chembl_id_lookup a, activities b, chembl28_ncats.spectral_hk c, assays d, compound_structures e
where a.entity_id = b.molregno
and b.molregno = e.molregno
and a.chembl_id = c.chembl_id
#and c.h1='WU7MZ9LTZ'
#and c.h2='NN51MHRR4DXYNAQ6H78'
#and c.h2='VXL4K9UW29154K6U6NH'
and c.h2='3U51T8HJB72X2TJJJC1'
and standard_relation='='
and b.assay_id = d.assay_id
and d.assay_id = 774053
;

select d.assay_id,d.description,type,units,count(distinct b.molregno) as N, avg(b.value), std(b.value)
from chembl_id_lookup a, activities b, chembl28_ncats.spectral_hk c, assays d
 where a.entity_id = b.molregno
and a.chembl_id = c.chembl_id
#and c.h1='WU7MZ9LTZ'
#and c.h2='NN51MHRR4DXYNAQ6H78'
#and c.h2='VXL4K9UW29154K6U6NH'
and c.h2='3U51T8HJB72X2TJJJC1'
and standard_relation='='
and b.assay_id = d.assay_id
group by d.assay_id,d.description,type,units
having count(*) > 1
order by N desc
;

select e.h2 as spectral_h2, d.pref_name as target, h.pref_name as target_class, count(*) as cnt, 
avg(-log(10,1e-9*b.standard_value)) as avg_act, std(-log(10,1e-9*b.standard_value)) as std_act
from chembl_id_lookup a, activities b, assays c, target_dictionary d, chembl28_ncats.spectral_hk e,
target_components f, component_class g, protein_classification h
where a.entity_id = b.molregno
and c.assay_id = b.assay_id
and c.tid = d.tid
and c.tid = f.tid
and f.component_id = g.component_id
and g.protein_class_id = h.protein_class_id
and a.chembl_id = e.chembl_id
and a.entity_type='COMPOUND'
and c.confidence_score >= 6
and standard_relation='='
and standard_type = 'Ki'
and d.target_type = 'SINGLE PROTEIN'
and d.tax_id=9606
group by e.h2, d.pref_name
having count(a.chembl_id) > 9
order by std_act, avg_act desc, e.h2, cnt desc
;
