function [surface1_R,surface2_R,thickness,surface1_ind,ya,ya_,yb,yb_,ua,ua_,ub,ub_,asphere,NA]=basic_data_generator(surface_num)
cvon% build up the connection
%cmd='res "C:\Users\32141\Downloads\CodeVv2007a\CodeVv2007a\test.len"';
%cvcmd(cmd)
%surface_num=2;
cvin;%import a .seq or a .len file from the cd of matlab
cmd = 'buf del b0; buf; THO; buf n';
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(0);%return the result as a double
line=1;
index=1;
if surface_num>=10
    index=2;
end
while 1
    if result(line,1:index)==num2str(surface_num)
        break
    end
    line=line+1;
end
surface_1=line;
%%%%%%%%aberration graber
database_temp=str2num(result(surface_1,:));
S11=database_temp(2);
S21=database_temp(3);
S31=database_temp(4);
S41=database_temp(5);
%%%%%%%%%%%%%%%%%%%%%
while 1
    if result(line,1:index)==num2str(surface_num+1)
        break
    end
    line=line+1;
end
surface_2=line;
database_temp=str2num(result(surface_2,:));
S12=database_temp(2);
S22=database_temp(3);
S32=database_temp(4);
S42=database_temp(5);
S1=S11+S12;
S2=S21+S22;
S3=S31+S32;
S4=S41+S42;
%%%%%%%%%%%%%%%%%%CURVATURE
cmd = sprintf("buf del b1; buf put b1 eva(rdy s%i)",surface_num);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
surface1_R=str2num(result(4:end));

cmd = sprintf("buf del b1; buf put b1 eva(rdy s%i)",surface_num+1);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
surface2_R=str2num(result(4:end));
%%%%%%%%%%%%%%%%%%%%%%%%index
cmd = sprintf("buf del b1; buf put b1 eva(ind s%i)",surface_num);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
surface1_ind=str2num(result(4:end));

cmd = sprintf("buf del b1; buf put b1 eva(ind s%i)",surface_num+1);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
surface2_ind=str2num(result(4:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THICKNESS
cmd = sprintf("buf del b1; buf put b1 eva(thi s%i)",surface_num);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
thickness=str2num(result(4:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%asphere coefficient
cmd = sprintf("buf del b1; buf put b1 eva((A s%i)*100000000)",surface_num);
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
asphere=str2num(result(4:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%NA
cmd = sprintf("buf del b1; buf put b1 eva(NA)");
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(1);
NA=str2num(result(4:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%angle, height
cmd="buf del b0;buf;FIO SO..I;buf n";
cvcmd(cmd);
result =cvbuf(0);
line=1;
while 1
    if result(line,1:index)==num2str(surface_num)
        break
    end
    line=line+1;
end
surface_1=line;
database_temp=str2num(result(surface_1,:));
ya=database_temp(2);

yb=database_temp(5);
%%%%%%%%%%%%u is before surface
database_temp=str2num(result(surface_1-1,:));
ua=database_temp(3);
ub=database_temp(6);

line=1;
while 1
    if result(line,1:index)==num2str(surface_num+1)
        break
    end
    line=line+1;
end
surface_2=line;
database_temp=str2num(result(surface_2,:));
ya_=database_temp(2);
yb_=database_temp(5);
%database_temp=str2num(result(surface_2-1,:));
ua_=database_temp(3);
ub_=database_temp(6);




cvoff % kills the connection