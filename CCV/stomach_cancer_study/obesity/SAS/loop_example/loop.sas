/* Calls to macros begin with % symbol */
/* Calls to macro variables begin with & symbol */
/* -------------------------------------------------------------------------- */
/* CODE BEGINS HERE */

/* Run in analyses2006/analise.sas under stratcod parameter, e.g.: */

stratcod =

    /* Assign list of variables to macro variable  */
    /* -------------------------------------------------------------------------- */
    %let input_var_list = var1 var2; /*list of variables, var3 var4 etc. ... */

    /* Create macro to run desired function through loop */
    /* -------------------------------------------------------------------------- */
    /* e.g. get boxplot for desired list of variables */

    %macro boxplot_loop(var_list = ,); 
        
        /* Check number of words (variables) in &var_list */
        %let len_var_list = %eval(%sysfunc(countw(&var_list.)));
        
        /* Iterate through sequence of variables */
        %do _i = 1 %to &len_var_list;

            /* Select i'th target variable and assign to &target_var */
            %let target_var = %scan(&var_list., &_i);

            /* Run desired code to loop, with reference to &target_var */
            /* e.g. produce png of boxplot of target variable by cohort 
            and name png file by target variable */
            ods graphics / imagename="Boxplot-&target_var" imagefmt=png;
            proc sgplot; 
                title "Boxplot-&target_var";
                vbox &target_var / category=dataname; * x-axis 
            run;
        %end;
    %mend boxplot_loop;

    /* Run macro function on list of all input variables  */
    %boxplot_loop(var_list = &input_var_list);

/* Alternatively, if you plan to reuse the boxplot loop for multiple scripts:
you can include the above 'boxplot_loop' macro in another file somewhere in 
your own directory, e.g. loop_ref.sas, then add '%include "loop_ref.sas";' at 
the top of your work script where you call %boxplot_loop . */
