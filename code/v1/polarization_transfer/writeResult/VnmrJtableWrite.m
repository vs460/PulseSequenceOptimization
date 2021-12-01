% function to write the result of INEPT optimization to real-time table files
% for VnmrJ
function[] = VnmrJtableWrite(protonAmp,protonPhase,protonTaus,carbonAmp,carbonPhase,carbonTaus)

% proton amplitude
fileID = fopen('vs_INEPT_protonAmps','w');
fprintf(fileID,'t1 = \n');
fprintf(fileID,'%3.1f, ',protonAmp);
fclose(fileID);
% proton phase
fileID = fopen('vs_INEPT_protonPhases','w');
fprintf(fileID,'t2 = \n');
fprintf(fileID,'%5.2f, ,',protonPhase);
fclose(fileID);
% proton taus
fileID = fopen('vs_INEPT_protonTaus','w');
fprintf(fileID,'t3 = \n');
fprintf(fileID,'%f, ',protonTaus);
fclose(fileID);
% carbon amplitude
fileID = fopen('vs_INEPT_carbonAmps','w');
fprintf(fileID,'t4 = \n');
fprintf(fileID,'%3.1f, ',carbonAmp);
fclose(fileID);
% carbon phase
fileID = fopen('vs_INEPT_carbonPhases','w');
fprintf(fileID,'t5 = \n');
fprintf(fileID,'%5.2f, ',carbonPhase);
fclose(fileID);
% carbon taus
fileID = fopen('vs_INEPT_carbonTaus','w');
fprintf(fileID,'t6 = \n');
fprintf(fileID,'%f, ',carbonTaus);
fclose(fileID);

end

