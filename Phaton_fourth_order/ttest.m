cvon% build up the connection
%cmd='res "C:\Users\32141\Downloads\CodeVv2007a\CodeVv2007a\test.len"';
%cvcmd(cmd)
surface_num=2;
cvin('test.len');%import a .seq or a .len file from the cd of matlab
cmd = 'buf del b0; buf; THO; buf n';
cvcmd(cmd);%similar to what you type in the command window
result =cvbuf(0);%return the result as a double
line=1;
while 1
    if result(line,1)==num2str(surface_num)
        break
    end
    line=line+1;
end
surface_1=line;
%%%%%%%%aberration graber
database_temp=str2num(result(surface_1,:))
S11=database_temp(2);
S21=database_temp(3);
S31=database_temp(4);
S41=database_temp(5);
%%%%%%%%%%%%%%%%%%%%%
while 1
    if result(line,1)==num2str(surface_num+1)
        break
    end
    line=line+1;
end
surface_2=line;
database_temp=str2num(result(surface_2,:))
S12=database_temp(2);
S22=database_temp(3);
S32=database_temp(4);
S42=database_temp(5);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%angle, height
cmd="buf del b0;buf;FIO SO..I;buf n";
cvcmd(cmd);
result =cvbuf(0);
line=1;
while 1
    if result(line,1)==num2str(surface_num)
        break
    end
    line=line+1;
end
surface_1=line;
database_temp=str2num(result(surface_1,:));
ya=database_temp(2);
ua=database_temp(3);

yb=database_temp(5);
ub=database_temp(6);

line=1;
while 1
    if result(line,1)==num2str(surface_num+1)
        break
    end
    line=line+1;
end
surface_2=line;
database_temp=str2num(result(surface_2,:));
ya_=database_temp(2);
ua_=database_temp(3);

yb_=database_temp(5);
ub_=database_temp(6);




cvoff % kills the connection