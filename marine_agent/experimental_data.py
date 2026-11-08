"""Strict, bounded experimental table ingestion. No model-generated measurements."""
from __future__ import annotations
import base64, csv, hashlib, io, json, math, statistics, zipfile
from collections import defaultdict

COLUMNS = ['candidate_key','metric','value','unit','qualifier','data_kind','sd','n','replicate_id',
           'substrate','temperature_c','ph','protocol','control','batch','construct','source']
# canonical unit and multiplier; no conversion between activity and catalytic constants.
UNITS = {
 'kcat': {'s^-1':('s^-1',1), 'min^-1':('s^-1',1/60)},
 'km': {'M':('mM',1000), 'mM':('mM',1), 'uM':('mM',.001), 'µM':('mM',.001)},
 'kcat_km': {'M^-1 s^-1':('M^-1 s^-1',1), 'M^-1 min^-1':('M^-1 s^-1',1/60)},
 'functional_tm': {'C':('°C',1),'°C':('°C',1)},
 'specific_activity': {'U/mg':('U/mg',1)}, 'relative_activity': {'%':('%',1)},
}
MISSING = {'','NA','N/A','null','None'}

def number(value, optional=False):
    if str(value).strip() in MISSING:
        if optional:return None
        raise ValueError('缺少数值')
    result=float(value)
    if not math.isfinite(result):raise ValueError('数值必须有限')
    return result

def read_table(filename, encoded):
    try: raw=base64.b64decode(encoded,validate=True)
    except Exception:raise ValueError('文件编码无效')
    if not raw or len(raw)>4_000_000:raise ValueError('文件须小于 4 MB')
    if filename.lower().endswith('.csv'):
        try: rows=list(csv.reader(io.StringIO(raw.decode('utf-8-sig')),strict=True))
        except Exception:raise ValueError('CSV 须为 UTF-8 编码')
    elif filename.lower().endswith('.xlsx'):
        from openpyxl import load_workbook
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                if len(z.infolist())>1000 or sum(i.file_size for i in z.infolist())>20_000_000:raise ValueError()
                if any('externallink' in i.filename.lower() for i in z.infolist()):raise ValueError()
            book=load_workbook(io.BytesIO(raw),read_only=True,data_only=False,keep_links=False)
            try:
                if len(book.worksheets)!=1:raise ValueError('仅支持一个工作表')
                rows=[]
                for row in book.active.iter_rows():
                    if len(rows)>10000 or len(row)>len(COLUMNS):raise ValueError('表格超过行列限制')
                    if any(c.data_type=='f' for c in row):raise ValueError('请将公式转换为值后再导入')
                    rows.append(['' if c.value is None else str(c.value) for c in row])
            finally:book.close()
        except ValueError:raise
        except Exception:raise ValueError('无法读取 Excel 文件')
    else:raise ValueError('仅支持 CSV 或 XLSX')
    rows=[r for r in rows if any(str(v).strip() for v in r)]
    if len(rows)<2 or len(rows)>10001:raise ValueError('需要 1–10000 条数据')
    header=[str(v).strip() for v in rows[0]]
    if len(header)!=len(set(header)) or set(header)!=set(COLUMNS):raise ValueError('列名不匹配，请使用模板')
    return raw,[dict(zip(header,r)) if len(r)==len(header) else {'_invalid':'列数不一致'} for r in rows[1:]]

def normalize_row(row, candidates):
    if '_invalid' in row:raise ValueError(row['_invalid'])
    row={k:str(v).strip() for k,v in row.items()}
    if any(len(v)>2000 for v in row.values()):raise ValueError('单元格内容过长')
    key=row['candidate_key']
    if key not in candidates:raise ValueError('候选编号不存在：'+key)
    metric=row['metric'];unit=row['unit']
    if metric not in UNITS or unit not in UNITS[metric]:raise ValueError('不支持的指标或单位：'+metric+' / '+unit)
    canonical,factor=UNITS[metric][unit]
    value=number(row['value'],True);sd=number(row['sd'],True);n=number(row['n'],True)
    qualifier=row['qualifier'] or ('missing' if value is None else '=')
    if qualifier not in {'=','<','>','missing'} or ((value is None)!=(qualifier=='missing')):raise ValueError('数值与限定符不一致')
    if value is not None and metric!='functional_tm' and value<0:raise ValueError('该指标不能为负数')
    if sd is not None and (sd<0 or qualifier!='='):raise ValueError('标准差只能用于精确测量值且不能为负')
    if n is not None and (n<1 or n!=int(n)):raise ValueError('n 须为正整数')
    kind=row['data_kind']
    if kind not in {'summary','replicate'}:raise ValueError('data_kind 须为 summary 或 replicate')
    if kind=='replicate' and (not row['replicate_id'] or sd is not None or n not in {None,1}):raise ValueError('重复测量须有 replicate_id，不能携带汇总标准差或 n>1')
    if kind=='summary' and row['replicate_id']:raise ValueError('汇总值不能冒充原始重复测量')
    if not row['source']:raise ValueError('必须填写实验来源')
    temp=number(row['temperature_c'],True);ph=number(row['ph'],True)
    if ph is not None and not 0<=ph<=14:raise ValueError('pH 超出 0–14')
    if temp is not None and not -273.15<temp<1000:raise ValueError('温度超出有效范围')
    conditions={k:(None if row[k] in MISSING or row[k].lower() in {'unknown','unreported','未知','未提供','未报告'} else row[k]) for k in ['substrate','protocol','control','batch','construct','source']}
    conditions.update(temperature_c=temp,ph=ph)
    return {'candidate_key':key,'candidate_id':candidates[key]['id'],'sequence_sha256':candidates[key]['sequence_sha256'],
            'metric':metric,'raw_value':row['value'],'raw_unit':unit,'value':value*factor if value is not None else None,
            'unit':canonical,'qualifier':qualifier,'data_kind':kind,'sd':sd*factor if sd is not None else None,
            'n':int(n) if n is not None else None,'replicate_id':row['replicate_id'] or None,'conditions':conditions,
            'provenance':'user_imported_experimental_table_not_independently_verified'}

def row_hash(row):return hashlib.sha256(json.dumps(row,sort_keys=True,ensure_ascii=False).encode()).hexdigest()

def comparisons(rows):
    groups=defaultdict(list)
    for row in rows:
        key=json.dumps([row['metric'],row['unit'],row['data_kind'],row['conditions']],sort_keys=True)
        groups[key].append(row)
    result=[]
    for values in groups.values():
        first=values[0];conditions=first['conditions'];issues=[]
        required=['substrate','ph','protocol','control','batch','construct']
        # Functional melting point is the fitted outcome of a temperature ramp.
        if first['metric']!='functional_tm':required.append('temperature_c')
        missing=[k for k in required if conditions.get(k) is None]
        if missing:issues.append('条件缺失：'+', '.join(missing))
        by_candidate=defaultdict(list)
        for row in values:by_candidate[row['candidate_key']].append(row)
        summaries=[]
        for key,measurements in by_candidate.items():
            exact=[r['value'] for r in measurements if r['qualifier']=='=']
            if len(exact)!=len(measurements):issues.append(key+' 有缺失或检出限数据，未纳入精确排序')
            if first['data_kind']=='summary' and len(measurements)>1:issues.append(key+' 存在多个汇总值，需要人工指定分析方案')
            if exact and len(exact)==len(measurements):
                summaries.append({'candidate_key':key,'value':statistics.mean(exact),
                                  'sd':statistics.stdev(exact) if first['data_kind']=='replicate' and len(exact)>1 else measurements[0]['sd'],
                                  'n':len(exact) if first['data_kind']=='replicate' else measurements[0]['n']})
        eligible=not issues and len(summaries)>1
        summaries.sort(key=(lambda r:r['value']) if eligible else (lambda r:r['candidate_key']))
        result.append({'metric':first['metric'],'unit':first['unit'],'data_kind':first['data_kind'],'conditions':conditions,
                       'comparable':eligible,'issues':list(dict.fromkeys(issues)), 'values':summaries,
                       'measurement_counts':{'exact':sum(r['qualifier']=='=' for r in values),'missing':sum(r['qualifier']=='missing' for r in values),'censored':sum(r['qualifier'] in {'<','>'} for r in values)},
                       'interpretation':'条件完整时仅按数值升序展示；高低是否有利取决于研究目标，不代表显著差异或工业表现。'})
    return result
