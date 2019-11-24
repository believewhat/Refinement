load('alarm');data=data';
LGObj = ConstructLGObj(data);
u = 4;
mex('XMY0.cpp');
[score,Order]=XMY0(LGObj,dag);
ns=LGObj.VarRangeLength;
DAG2=learn_struct_K2(data',ns,Order,'scoring_fn','bic');
h = view(biograph( DAG2 ));
h=view(biograph( dag ));
cnt=0;
for i=1:37
    for j=i+1:37
        if dag(i,j)~=DAG2(i,j) cnt=cnt+1;end
    end
end
cnt
DAG3=learn_struct_K2(data',ns,randperm(37),'scoring_fn','bic');
cnt=0;
for i=1:37
    for j=i+1:37
        if DAG3(i,j)~=dag(i,j)==0
            cnt=cnt+1;
        end
    end
end
cnt
m=size(anss,2);
vis=zeros(1,20000);
for i=1:m
    vis(anss(i)+1)=1;
end
for i=1:m
    datt(i,:)=data(1:6,anss(i)+1)';
    trainlabel(i)=2^(data(7,anss(i)+1)-1)+data(8,anss(i)+1)-1;
end
kk=0;
while(kk<5000)
    i=ceil(rand*20000);
    if vis(i)==1
        continue;
    end
    vis(i)=1;
    testdata(int32(kk+1),:)=data(1:6,i)';
    testlabel(int32(kk+1))=2^(data(7,i)-1)+data(8,i)-1;
    kk=kk+1;
end
nb =ClassificationKNN.fit(datt, trainlabel,'NumNeighbors',3);
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
for i=1:20000
    datt(i,:)=data(1:6,i)';
    trainlabel(i)=2^(data(7,i)-1)+data(8,i)-1;
end
clear datt trainlabel accuracy nb predict_label nb predict_label 
nb =svmtrain(datt, trainlabel);
predict_label=svmclassify(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;

ns=LGObj.VarRangeLength;
DAG2=learn_struct_K2(datt',ns,order,'scoring_fn','bic');
h = view(biograph( dag ));
h2=view(biograph( DAG2 ));
data=xlsread('10.xlsx');
data=data+1;
mex('k2change.cpp');
dag=zeros(91,91);
[DAG3,Order]=k2change(LGObj,dag);
ns=LGObj.VarRangeLength;
load('asia');
LGObj = ConstructLGObj(data');
u = 4;
mex('XMY.cpp');
[score,Order]=XMY(LGObj,dag);
DAG2=learn_struct_K2(data,ns,Order,'scoring_fn','bic');
h = view(biograph( DAG2 ));
cnt=0;
for i=1:37
    for j=1:37
        if DAG2(i,j)>0 && dag(i,j)==0 && dag(j,i)==1
            cnt=cnt+1;
        end
    end
end
cnt
h = view(biograph( DAG2 ));

h = view(biograph( dag ));
%correct=cal(8,168,9,600,9,DAG2,data(601:768,1:9),data(1:600,1:9))
bnet2 = bayes_update_params(bnet,data);

x=[0.5 0.7 0.2 0.4 2.5 1.5 -0.2 -0.5 0.1];
y=[1 2 1 1 1 1 1 1 1];








load('alarm');
LGObj = ConstructLGObj(data(1:36,1:20000)');
u = 4;
mex('XMY2.cpp');
[score,anss]=XMY2(LGObj);
m=size(anss,2);
vis=zeros(1,20000);
for i=1:m
    vis(anss(i)+1)=1;
end
for i=1:m
    datt(i,:)=data(:,anss(i)+1)';
    trainlabel(i)=data(37,anss(i)+1);
end
kk=0;
while(kk<5000)
    i=ceil(rand*20000);
    if vis(i)==1
        continue;
    end
    vis(i)=1;
    testdata(int32(kk+1),:)=data(1:36,i)';
    testlabel(int32(kk+1))=data(37,i);
    kk=kk+1;
end
yfit=trainclass.predictFcn(testdata);
nb =ClassificationKNN.fit(datt, trainlabel,'NumNeighbors',6);
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
accuracy=length(find(yfit'==testlabel))/length(testlabel)*100;
kk=1;
for i=1:20000
    datt(kk,:)=data(1:36,i)';
    trainlabel(kk)=data(37,i);
    kk=kk+1;
end
clear datt trainlabel accuracy nb predict_label nb predict_label


[model_pos,model_neg] = FindGuassianModel(datt,trainlabel);
for i=1:1000
    a(i)=chi2inv(0.95,i);
end
a=a';
save kafa.txt a -ascii -double