LGObj = ConstructLGObj(data');
u = 4;
mex('XMY.cpp');
[score,Order]=XMY(LGObj,dag);
DAG2=learn_struct_K2(LGObj.VarSample,ns,Order,'scoring_fn','bic');
h = view(biograph( DAG2 ));