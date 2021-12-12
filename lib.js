
const prefixSums = vs => vs.reduce((acc,v,i)=>{acc[i]=(i?acc[i-1]:0)+v;return acc},[])
const sum = vs => vs.reduce((a,b)=>a+b,0)
const avg = vs => sum(vs)/vs.length;
const max = vs => Math.max.apply(Math,vs);
const min = vs => Math.min.apply(Math,vs);
// variance, d2
const avg_err2 = (n, sum, sum2) => sum2/n - (sum/n)**2
const entropy = (n,sum,sum2) => n*0.5*(1 + Math.log(avg_err2(n,sum,sum2)*2*Math.PI))
const plus = (a,b) => ({x:a.x+b.x,y:a.y+b.y});
const minus = (a,b) => ({x:a.x-b.x,y:a.y-b.y});
const scale = (s,p) => ({x:s*p.x,y:s*p.y});
const scalar = (a,b) => a.x*b.x + a.y*b.y;
const always = (x) => () => x;
const length = p => Math.sqrt(p.x**2 + p.y**2)
const distance = (a,b) => length(minus(a,b));
const normalize = p => scale(1.0/length(p),p);
const deg2rad = d => Math.PI*d/180;
const polar = (a,r) => ({x:Math.cos(a)*r,y:Math.sin(a)*r});
const flatten = ass => [].concat.apply([],ass);
function split(vs){
  const sorted=vs.slice(0).sort((a,b)=>a-b);
  const sum = prefixSums(sorted);
  const sum2 = prefixSums(vs.map((v)=>v*v));
  let bestCost = +Infinity;
  let best;
  let bestI;
  sorted.forEach((v,i)=>{
    if(i&&v!=sorted[i-1]){
      const cost = entropy(i,sum[i-1],sum2[i-1])+entropy(sorted.length-i,sum.at(-1)-sum[i-1],sum2.at(-1)-sum2[i-1]);
      if(cost<bestCost){
        bestCost=cost;
        best=(v+sorted[i-1])/2;
        bestI=i;
      }
    }
  })
  const res={threshold:best, low: sum[bestI-1]/bestI, hi: (sum.at(-1)-sum[bestI-1])/(sorted.length-bestI),bestCost};
  res.mid = (res.low + res.hi)/2;
  return res;
}
function modDist(a,b,m){
  const d=Math.abs(a-b);
  return Math.min(d,m-d);
}
function optimalLine(points){
  /*
  minimize: sum(p:points){(p.x*cos(a) + p.y*sin(a) - off)^2}
  minimize: sum(p:points){(f(a,p) - avg(p:points){f(a,p)})^2}
  minimize: variance(p:points){f(a,p)}
  minimize: sum(p:points){f(a,p)^2}/n - sum(p:points){f(a,p)}^2/n 
  minimize: sum(p:points){(p.x*cos(a) + p.y*sin(a))^2}/n - sum(p:points){p.x*cos(a) + p.y*sin(a)}^2/n^2 
  minimize: sum(p:points){p.x^2*cos(a)^2 + p.y^2*sin(a)^2 + 2*p.x*p.y*sin(a)*cos(a)}/n - (sum(p:points){p.x*cos(a)}+sum(p:points){p.y*sin(a)})^2/n^2
  minimize: sum(p:points){p.x^2*cos(a)^2 + p.y^2*sin(a)^2 + 2*p.x*p.y*sin(a)*cos(a)}/n - (cos(a)sum(p:points){p.x}+sin(a)sum(p:points){p.y})^2/n^2 
  minimize: Sx2/n*cos(a)^2 + Sy2/n*sin(a)^2 + 2*Sxy/n*sin(a)*cos(a) - (cos(a)Sx+sin(a)Sy)^2/n^2 
  minimize: Sx2/n*cos(a)^2 + Sy2/n*sin(a)^2 + 2*Sxy/n*sin(a)*cos(a) - (Sx^2*cos(a)^2+Sy^2sin(a)^2+2*Sx*Sy*cos(a)sin(a))/n^2 
    t=cos(a)^2
    sin(a)^2=(1-t)
    cos(a)*sin(a)= +- sqrt(t(1-t))
  minimize: Sx2/n*t + Sy2/n*(1-t) + 2*Sxy/n*+-sqrt(t(1-t)) - (Sx^2*t+Sy^2(1-t)+2*Sx*Sy*+-sqrt(t(1-t)))/n^2 


  diff by t of: sqrt(t(1-t)) 
  diff by t of: (t(1-t))^0.5
            is: 0.5*(t(1-t))^-0.5 * (t(1-t))'
            is: 0.5*(t(1-t))^-0.5 * ((1-t)-t)
            is: 0.5*(t(1-t))^-0.5 * (1-2t)

  diff by t of: Sx2/n*t + Sy2/n*(1-t) + 2*Sxy/n*+-sqrt(t(1-t)) - (Sx^2*t+Sy^2(1-t)+2*Sx*Sy*+-sqrt(t(1-t)))/n^2 
            is: Sx2/n   - Sy2/n       + 2*Sxy/n*+-0.5*(t(1-t))^-0.5 * (1-2t) - (Sx^2 - Sy^2 + 2*Sx*Sy*+-0.5*(t(1-t))^-0.5 * (1-2t))/n^2 
            is: Sx2/n   - Sy2/n       + +-Sxy/n*(t(1-t))^-0.5 * (1-2t) - (Sx^2 - Sy^2 + +-Sx*Sy*(t(1-t))^-0.5 * (1-2t))/n^2 
            is: Sx2/n   - Sy2/n       - (Sx^2 - Sy^2)/n^2   + +-(Sxy/n - Sx*Sy/n^2)*(t(1-t))^-0.5 * (1-2t)  
  should be zero, so:
  Sx2/n   - Sy2/n       - (Sx^2 - Sy^2)/n^2   =  -+ (Sxy/n - Sx*Sy/n^2)*(t(1-t))^-0.5 * (1-2t)  
  (Sx2/n - Sy2/n - (Sx^2 - Sy^2)/n^2) / (Sxy/n - Sx*Sy/n^2)   =  +- (t(1-t))^-0.5 * (2t-1)  
  C   =  +- (t(1-t))^-0.5 * (2t-1)  
  C**2 * t(1-t)  =   (2t-1)**2  

      -C**2 * t**2   + C**2 * t =    4 * t**2 - 4 * t + 1

    (-4-C**2) * t**2 + (C**2 + 4) * t - 1  = 0;
     t**2 -  t + 1/(C**2 + 4)  = 0;
  */
  const Sx=sum(points.map(p=>p.x))
  const Sy=sum(points.map(p=>p.y))
  const Sx2=sum(points.map(p=>p.x**2))
  const Sy2=sum(points.map(p=>p.y**2))
  const Sxy=sum(points.map(p=>p.x*p.y))
  const variance = a =>{
    const c=Math.cos(a);
    const s=Math.sin(a);
    return Sx2/n*c**2 + Sy2/n*s**2 + 2*Sxy/n*s*c - (Sx**2*c**2+Sy**2*s**2+2*Sx*Sy*c*s)/n**2;
  }

  const n = points.length;
  const C = ((Sx2 - Sy2)/n - (Sx*Sx-Sy*Sy)/n**2)/(Sxy/n - Sx*Sy/n**2);
  const a= 1;
  const b = -1;
  const c = 1/(C**2+4);
  const delta = b**2 - 4*a*c; // = 1-4/(C**2+4), should be 0 for C=0 and larger than 0 for any other C
  console.assert(0<=delta);
  const ts = [1,-1].map(sign => (-b + sign*Math.sqrt(delta))/(2*a));
  // so, we have two different ts,
  // and each t is cos(a)^2, so maps to two different cos(a)
  // finally each cos(a) maps to two different angles in range [0,2PI), but we want something in range [0,PI), and in this range it's unique
  const candidate_angles=flatten(ts.filter(t=> 0<=t && t<=1).map(t => [1,-1].map(sign => Math.acos(sign*Math.sqrt(t)))));
  const ang = candidate_angles.map(a => [variance(a),a]).sort((a,b)=>a[0]-b[0])[0][1];
  const off = (Sx*Math.cos(ang)+Sy*Math.sin(ang))/n;
  return {a:ang,off,x:Math.cos(ang),y:Math.sin(ang)};
}
function solve(A,b){
  A=A.map(row=>row.slice(0));
  b=b.slice(0)
  //Ax=b
  const vars_num=A[0].length;
  const eq_num=A.length;
  for(let c=0;c<vars_num;++c){
    let best_r = c;
    for(let r=c+1;r<eq_num;++r){
      if(Math.abs(A[best_r][c])<Math.abs(A[r][c])){
        best_r=r;
      }
    }
    {
      let r=c;
      let tmp_row=A[r];
      A[r]=A[best_r];
      A[best_r]=tmp_row;
      let tmp=b[r];
      b[r]=b[best_r];
      b[best_r]=tmp;
    }
    for(let r=0;r<eq_num;++r)if(r!=c){
      const s=A[r][c]/A[c][c];
      for(let i=0;i<vars_num;++i){
        A[r][i]-= A[c][i]*s;
      }
      b[r]-=b[c]*s;
    }
  }
  let solution=[];
  for(let v=0;v<vars_num;++v){
    solution[v]=b[v]/A[v][v]
  }
  return solution;
}
function evaluate(a,x){
  if(!Number.isFinite(x)){
    return a.at(-1)*(x**(a.length-1));
  }
  return a.reduce((acc,v,i)=>acc+v*(x**i),0);
}
function zign(x){
  return x<-1e-9?-1:1e-9<x?1:0;
}
function bisect_zero(a,min,max){
  const increasing = evaluate(a,min)<=evaluate(a,max);
  const larger = (x,dir) => (Math.abs(x)*2+1)*Math.sign(dir);
  const split = (min,max)=>Number.isFinite(min)?
      (Number.isFinite(max) ? (max+min)*0.5: larger(min,max)):
      (Number.isFinite(max) ? larger(max,min) : 0);
  for(let i=0;i<10000;++i){
    const mid=split(min,max);
    const val=evaluate(a,mid);
    if(!zign(val)){
      return mid;
    }
    if((val<=0) == increasing ){
      min=mid;
    }else{
      max=mid;
    }
    if(i==10000-1){
      debugger;
    }
  }
  return min;
}
function derivative(a){
  return a.slice(1).map((v,k)=>(k+1)*v)
}
function zeros(a){
  while(a.length && !a.at(-1))a.pop()
  if(!a.length){
    return [0];
  }
  if(a.length==1){
    return a[0]?[]:[0];
  }
  const found=[];
  [-Infinity,...zeros(derivative(a)),+Infinity].forEach((x,i,xs)=>{
    if(i&&zign(evaluate(a,x))!=zign(evaluate(a,xs[i-1]))){
      found.push(bisect_zero(a,xs[i-1],x));
    }
  })
  return found;
}
function boxesIntesect(a,b){
  return a.min.x <= b.max.x && b.min.x <= a.max.x &&
         a.min.y <= b.max.y && b.min.y <= a.max.y;
}