[num,raw] = xlsread('alarm.xlsx');
load('alarm');
[n,m]=size(raw);
for i=2:m
    cnt=1;
    map=containers.Map();
    for j=2:n
        if isKey(map,char(raw(j,i)))==0
            map(char(raw(j,i)))=cnt;
            data(j-1,i-1)=cnt;
            cnt=cnt+1;
        else
            data(j-1,i-1)=map(char(raw(j,i)));
        end
    end
end

LGObj = ConstructLGObj(data);
mex('XMY2.cpp');
[score,anss]=XMY2(LGObj,dag,19980);
m=size(anss,2);
anss=anss+1;
vis=zeros(1,20000);
for i=1:m
    vis(anss(i))=1;
end
for i=1:m
    datt(i,:)=[data(anss(i),1:9) data(anss(i),11:37)];
    trainlabel(i)=data(anss(i),10);
end
kk=0;
while(kk<10000)
    i=ceil(rand*20000);
    if vis(i)==1
        continue;
    end
    vis(i)=1;
    testdata(int32(kk+1),:)=[data(i,1:9) data(i,11:37)];
    testlabel(int32(kk+1))=data(i,10);
    kk=kk+1;
end
nb =fitensemble(datt,trainlabel,'AdaBoostM2' ,500,'tree','type','classification');
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
nb =fitcecoc(datt,trainlabel);
t= templateSVM('Standardize',1);
nb =fitcecoc(datt,trainlabel,'Learners',t);
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
nb =fitcknn(datt, trainlabel);
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
nb =fitcnb(datt, trainlabel);
predict_label=nb.predict(testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
for i=1:20000
    datt(kk,:)=[data(i,1:9) data(i,11:37)];
    trainlabel(kk)=data(i,10);
    kk=kk+1;
end
vis=zeros(1,20000);
kk=0;
while(kk<50)
    i=floor(rand*20000);
    if vis(i)==1
        continue;
    end
    vis(i)=1;
    datt(int32(kk+1),:)=[data(i,1:9) data(i,11:37)];
    trainlabel(int32(kk+1))=data(i,10);
    kk=kk+1;
end
mex('calu.cpp');
[score]=calu(LGObj,dag,datt);
clear datt trainlabel accuracy nb predict_label nb predict_label nb predict_label

for i=1:36
i
var(testdata(:,i))
end



[num,raw] = xlsread('datawin95.csv');
[n,m]=size(raw);
for i=1:m
    cnt=1;
    map=containers.Map();
    for j=2:n
        if isKey(map,char(raw(j,i)))==0
            map(char(raw(j,i)))=cnt;
            data(j-1,i)=cnt;
            cnt=cnt+1;
        else
            data(j-1,i)=map(char(raw(j,i)));
        end
    end
end

LGObj = ConstructLGObj(data);
dag=zeros(76,76);
for i=1:76
    j=2;
    while j<=8 && ~isnan(nodeparnt(i,j))
        dag(nodeparnt(i,j)+1,i)=1;j=j+1;
    end
end
node=zeros(1,76);knode=1;
for i=1:76
    flag=0;
    for j=1:76
        if i==j continue;end
        if(dag(i,j)) flag=1;break;end
    end
    if ~flag node(knode)=i;knode=knode+1;end
end
%mex('XMY2.cpp');
maxdata=100;
[score,anss]=XMY2(LGObj,dag,20000-maxdata);
m=size(anss,2);
anss=anss+1;
vis=zeros(1,20000);
for i=1:m
    vis(anss(i))=1;
end
accuracy=100;accuracy3=100;accuracy4=100;accuracy5=100;
for h=1:knode-1
    tem1=node(h);
    for i=1:m
        datt(i,:)=[data(anss(i),1:tem1-1) data(anss(i),tem1+1:76)];
        trainlabel(i)=data(anss(i),tem1);
    end
    kk=0;
    vis2=vis;
    while(kk<10000)
        i=ceil(rand*20000);
        if vis2(i)==1
            continue;
        end
        vis2(i)=1;
        testdata(int32(kk+1),:)=[data(i,1:tem1-1) data(i,tem1+1:76)];
        testlabel(int32(kk+1))=data(i,tem1);
        kk=kk+1;
    end
    nb =fitensemble(datt,trainlabel,'AdaBoostM1' ,500,'tree','type','classification');
    predict_label=predict(nb,testdata);
    accuracy2=length(find(predict_label'==testlabel))/length(testlabel)*100;
    accuracy=min(accuracy,accuracy2);
    nb =fitcecoc(datt,trainlabel);
    t= templateSVM('Standardize',1);
    nb =fitcecoc(datt,trainlabel,'Learners',t);
    predict_label=predict(nb,testdata);
    accuracy2=length(find(predict_label'==testlabel))/length(testlabel)*100;
    accuracy4=min(accuracy4,accuracy2);
    clear datt trainlabel accuracy2 nb predict_label nb predict_label nb predict_label
    vis3=zeros(1,20000);
    kk=0;
    while(kk<maxdata)
        i=floor(rand*20000);
        if vis3(i)==1
            continue;
        end
        vis3(i)=1;
        datt(int32(kk+1),:)=[data(i,1:tem1-1) data(i,tem1+1:76)];
        trainlabel(int32(kk+1))=data(i,tem1);
        kk=kk+1;
    end
    kk=0;vis2=vis3;
    while(kk<10000)
        i=ceil(rand*20000);
        if vis2(i)==1
            continue;
        end
        vis2(i)=1;
        testdata(int32(kk+1),:)=[data(i,1:tem1-1) data(i,tem1+1:76)];
        testlabel(int32(kk+1))=data(i,tem1);
        kk=kk+1;
    end
    nb =fitensemble(datt,trainlabel,'AdaBoostM1' ,500,'tree','type','classification');
    predict_label=predict(nb,testdata);
    accuracy2=length(find(predict_label'==testlabel))/length(testlabel)*100;
    accuracy3=min(accuracy3,accuracy2);
    nb =fitcecoc(datt,trainlabel);
    t= templateSVM('Standardize',1);
    nb =fitcecoc(datt,trainlabel,'Learners',t);
    predict_label=predict(nb,testdata);
    accuracy2=length(find(predict_label'==testlabel))/length(testlabel)*100;
    accuracy5=min(accuracy5,accuracy2);
end
accuracy
accuracy3
accuracy4
accuracy5
vis=zeros(1,20000);
kk=0;
while(kk<50)
    i=floor(rand*20000);
    if vis(i)==1
        continue;
    end
    vis(i)=1;
    datt(int32(kk+1),:)=[data(i,1:29) data(i,31:76)];
    trainlabel(int32(kk+1))=data(i,30);
    kk=kk+1;
end
nb =fitcecoc(datt,trainlabel);
t= templateSVM('Standardize',1);
nb =fitcecoc(datt,trainlabel,'Learners',t);
predict_label=predict(nb,testdata);
accuracy=length(find(predict_label'==testlabel))/length(testlabel)*100;
clear datt trainlabel accuracy nb predict_label nb predict_label nb predict_label
