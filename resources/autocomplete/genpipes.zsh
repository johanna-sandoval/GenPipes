#compdef genpipes

# AUTOMATCALLY GENERATED by `shtab`


_shtab_genpipes_commands() {
  local _commands=(
    "ampliconseq:Version\: 4.1.2"
    "chipseq:Version\: 4.1.2"
    "covseq:Version\: 4.1.2"
    "dnaseq:Version\: 4.1.2"
    "dnaseq_high_coverage:Version\: 4.1.2"
    "epiqc:Version\: 4.1.2"
    "hicseq:Version\: 4.1.2"
    "methylseq:Version\: 4.1.2"
    "nanopore:Version\: 4.1.2"
    "nanopore_covseq:Version\: 4.1.2"
    "rnaseq:Version\: 4.1.2"
    "rnaseq_denovo_assembly:Version\: 4.1.2"
    "rnaseq_light:Version\: 4.1.2"
    "tumor_pair:Version\: 4.1.2"
  )
  _describe 'genpipes commands' _commands
}

_shtab_genpipes_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
)

_shtab_genpipes_ampliconseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[AmpliconSeq analysis type]:protocol:(qiime dada2)"
)

_shtab_genpipes_chipseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[Type of pipeline (default chipseq)]:protocol:(chipseq atacseq)"
)

_shtab_genpipes_covseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_dnaseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[DNAseq analysis type]:protocol:(mugqic mpileup light sv)"
)

_shtab_genpipes_dnaseq_high_coverage_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_epiqc_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_hicseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-e,--enzyme}"[Restriction Enzyme used to generate Hi-C library (default DpnII)]:enzyme:(DpnII HindIII NcoI MboI Arima)"
  {-t,--type}"[Hi-C experiment type (default hic)]:protocol:(hic capture)"
)

_shtab_genpipes_methylseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_nanopore_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_nanopore_covseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[Type of CoVSeQ analysis,basecalling on\/off (default without basecalling)]:protocol:(default basecalling)"
)

_shtab_genpipes_rnaseq_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[RNAseq analysis type]:protocol:(stringtie cufflinks)"
)

_shtab_genpipes_rnaseq_denovo_assembly_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-t,--type}"[RNAseq analysis type]:protocol:(trinity seq2fun)"
)

_shtab_genpipes_rnaseq_light_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
)

_shtab_genpipes_tumor_pair_options=(
  "(- :)"{-h,--help}"[show this help message and exit]"
  {-c,--config}"[config INI-style list of files\; config parameters are overwritten based on files order]:config:"
  {-s,--steps}"[step range e.g. \'1-5\', \'3,6,7\', \'2,4-8\']:steps:"
  {-o,--output-dir}"[output directory (default\: current)]:output_dir:"
  {-j,--job-scheduler}"[job scheduler type (default\: slurm)]:job_scheduler:(pbs batch daemon slurm)"
  {-f,--force}"[force creation of jobs even if up to date (default\: false)]"
  "--no-json[do not create JSON file per analysed sample to track the analysis status (default\: false i.e. JSON file will be created)]"
  "--report[create \'pandoc\' command to merge all job markdown report files in the given step range into HTML, if they exist\; if --report is set, --job-scheduler, --force, --clean options and job up-to-date status are ignored (default\: false)]"
  "--clean[create \'rm\' commands for all job removable files in the given step range, if they exist\; if --clean is set, --job-scheduler, --force options and job up-to-date status are ignored (default\: false)]"
  "--container[Run inside a container providing a validsingularity image path]:container:"
  {--genpipes_file,-g}"[Command file output path. This is the command used to process the data, or said otherwise, this command will \"run the Genpipes pipeline\". Will be redirected to stdout if the option is not provided.]:genpipes_file:"
  {-l,--log}"[log level (default\: info)]:log:(debug info warning error critical)"
  "--sanity-check[run the pipeline in \`sanity check mode\` to verify that all the input files needed for the pipeline to run are available on the system (default\: false)]"
  "--wrap[Path to the genpipe cvmfs wrapper script.
Default is genpipes\/ressources\/container\/bin\/container_wrapper.sh. This is a convenience options for using genpipes in a container]:wrap:"
  {-r,--readsets}"[readset file]:readsets_file:"
  {-d,--design}"[design file]:design_file:"
  "(- :)"{-v,--version}"[show the version information and exit]"
  {-p,--pairs}"[pairs file]:pairs:"
  "--profyle[adjust deliverables to PROFYLE folder conventions (Default\: False)]"
  {-t,--type}"[Tumor pair analysis type]:protocol:(fastpass ensemble sv)"
)


_shtab_genpipes() {
  local context state line curcontext="$curcontext"

  _arguments -C $_shtab_genpipes_options \
    ': :_shtab_genpipes_commands' \
    '*::: :->genpipes'

  case $state in
    genpipes)
      words=($line[1] "${words[@]}")
      (( CURRENT += 1 ))
      curcontext="${curcontext%:*:*}:_shtab_genpipes-$line[1]:"
      case $line[1] in
        ampliconseq) _arguments -C $_shtab_genpipes_ampliconseq_options ;;
        chipseq) _arguments -C $_shtab_genpipes_chipseq_options ;;
        covseq) _arguments -C $_shtab_genpipes_covseq_options ;;
        dnaseq) _arguments -C $_shtab_genpipes_dnaseq_options ;;
        dnaseq_high_coverage) _arguments -C $_shtab_genpipes_dnaseq_high_coverage_options ;;
        epiqc) _arguments -C $_shtab_genpipes_epiqc_options ;;
        hicseq) _arguments -C $_shtab_genpipes_hicseq_options ;;
        methylseq) _arguments -C $_shtab_genpipes_methylseq_options ;;
        nanopore) _arguments -C $_shtab_genpipes_nanopore_options ;;
        nanopore_covseq) _arguments -C $_shtab_genpipes_nanopore_covseq_options ;;
        rnaseq) _arguments -C $_shtab_genpipes_rnaseq_options ;;
        rnaseq_denovo_assembly) _arguments -C $_shtab_genpipes_rnaseq_denovo_assembly_options ;;
        rnaseq_light) _arguments -C $_shtab_genpipes_rnaseq_light_options ;;
        tumor_pair) _arguments -C $_shtab_genpipes_tumor_pair_options ;;
      esac
  esac
}



typeset -A opt_args
_shtab_genpipes "$@"
