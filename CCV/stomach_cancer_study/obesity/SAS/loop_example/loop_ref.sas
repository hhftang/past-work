filename sasmac "!SASROOT/sasautos";
options mautosource mrecall sasautos=(sasmac);

%macro boxplot_loop(var_list = ,); 
    %let len_var_list = %eval(%sysfunc(countw(&var_list.)));
    %do _i = 1 %to &len_var_list;
        %let target_var = %scan(&var_list., &_i);
        ods graphics / imagename="Boxplot-&target_var" imagefmt=png;
        proc sgplot; 
            title "Boxplot-&target_var";
            vbox &target_var / category=dataname;
        run;
    %end; 
%mend boxplot_loop;

/* -------- */
/* In your sas work script: */
/*

%include "~/loop_ref.sas";
...
stratcod =
    %let input_var_list = var1 var2; * etc.
    %boxplot_loop(var_list = &input_var_list);
,
...
*/