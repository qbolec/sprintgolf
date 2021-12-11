
const prefixSums = vs => vs.reduce((acc,v,i)=>{acc[i]=(i?acc[i-1]:0)+v;return acc},[])
const avg = vs => vs.reduce((a,b)=>a+b,0)/vs.length;
const max = vs => Math.max.apply(Math,vs);
const min = vs => Math.min.apply(Math,vs);
// variance, d2
const avg_err2 = (n, sum, sum2) => sum2/n - (sum/n)**2
const entropy = (n,sum,sum2) => n*0.5*(1 + Math.log(avg_err2(n,sum,sum2)*2*Math.PI))

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
