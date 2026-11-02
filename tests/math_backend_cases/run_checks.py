#!/usr/bin/env python3
"""运行包内四GL和两个原始GLML样本，检查有限输入与既有科学输出哈希。"""
from pathlib import Path
import argparse,json,subprocess,os,hashlib,tempfile
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--bin-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--only-gl',action='store_true');args=p.parse_args();binary=args.bin_dir.resolve();args.output.mkdir(parents=True,exist_ok=True)
manifest=json.loads((ROOT/'manifest.json').read_text())
for n,digest in manifest['input_sha256'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==digest,n
results=[]
for c in json.loads((ROOT/'cases.json').read_text()):
 if args.only_gl and c['target']!='aster_iod_gl':continue
 template=ROOT/c['setting'];lines=[]
 for line in template.read_text().splitlines():
  if '=' in line:
   k,v=map(str.strip,line.split('=',1))
   if k in ['file_observation','file_observer','file_EOP','file_initial_orbit']:line=k+' = '+str((template.parent/v).resolve())
  lines.append(line)
 d=Path(tempfile.mkdtemp(prefix=c['id']+'-',dir=args.output.resolve()));(d/'setting.txt').write_text('\n'.join(lines)+'\n');env=dict(os.environ,OMP_NUM_THREADS='4' if c['target']=='aster_iod_glml' else '1',OMP_DYNAMIC='FALSE')
 with (d/'console.txt').open('w') as stream:r=subprocess.run([str(binary/c['target']),'-I',str(d/'setting.txt')],cwd=d,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1200)
 assert r.returncode==c['exit'],(c['id'],r.returncode)
 for n,digest in c['references'].items():assert hashlib.sha256((d/n).read_bytes()).hexdigest()==digest,(c['id'],n)
 results.append({'case':c['id'],'expected_exit':r.returncode,'scientific_outputs_match':True,'directory':str(d)});print(json.dumps(results[-1]),flush=True)
(args.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
