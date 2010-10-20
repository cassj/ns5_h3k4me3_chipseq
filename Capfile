###
# Align the data with Bowtie and dump the resulting BAM files
# on S3
# Run Macs peak finding on them.

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, 'ami-20794c54'  #EC2 eu-west-1 64bit Lucid
set :instance_type, 'm1.large'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'

#note to self
#ami-52794c26 32bit Lucid
#ami-505c6924 64bit Maverick
#ami-20794c54 64bit Lucid

set :nhosts, 1
set :group_name, 'ns5_h3k4me3_chipseq'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 50  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1b'  #is where the ubuntu ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

set :input, "#{mount_point}/CMN038.export.txt"
set :ip, "#{mount_point}/CMN039.export.txt"

#grab the snapshot ID for the raw data (fastq)
task :get_snap_id, :roles=>:master do
  `curl http://github.com/cassj/ns5_h3k4me3_chipseq/raw/master/data/SNAPID > SNAPID `
end 


#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:mount_xfs


# Bowtie v0.12.7
bowtie_link = 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie%2F0.12.7%2F&ts=1287377285&use_mirror=mesh'

#should really save this as an ami, but no time at the moment
desc "install bowtie"
task :install_bowtie, :roles => group_name do
  run "sudo apt-get update"
  run "sudo apt-get install -y zip unzip"
  run "cd #{working_dir} && wget -Obowtie.zip #{bowtie_link}" 
  run "cd #{working_dir} && unzip bowtie.zip"
  run "sudo cp #{working_dir}/bowtie*/bowtie* /usr/local/bin/"
end 
before 'install_bowtie', 'EC2:start'




# s3 setup
desc "install s3 client"
task :install_s3, :roles => group_name do
  sudo 'apt-get update'
  sudo 'apt-get install -y s3cmd'
end
before 'install_s3', 'EC2:start'

desc "upload the s3 config file"
task :s3_config, :roles => group_name do
  upload(s3cfg, '/home/ubuntu/.s3cfg')
end 
before 's3_config', 'EC2:start'


#get the current mouse genome (which I already have on S3).
task :fetch_genome, :roles => group_name do
  run "s3cmd get --force s3://bowtie.mm9/mm9.ebwt.zip #{working_dir}/bowtie-0.12.7/indexes/mm9.ebwt.zip"
  run "rm -Rf #{working_dir}/bowtie-0.12.7/indexes/chr*"
  run "cd  #{working_dir}/bowtie-0.12.7/indexes && unzip -o mm9.ebwt.zip"
end 
before "fetch_genome","EC2:start"


# convert the export.txt file to fastq for alignment
task :make_fastq, :roles => group_name do 
  run "curl http://github.com/cassj/ns5_h3k4me3_chipseq/raw/master/scripts/export2fastq.pl > #{working_dir}/export2fastq.pl"
  run "chmod +x #{working_dir}/export2fastq.pl"
  run "sudo mv #{working_dir}/export2fastq.pl /usr/local/bin/"
  ip_o = ip.sub('.export.txt', '.fastq')
  input_o = input.sub('.export.txt', '.fastq')
  run "export2fastq.pl #{input} > #{input_o}"
  run "export2fastq.pl #{ip} > #{ip_o}"
end
before 'make_fastq', 'EC2:start' 



# run bowtie on the fastq file
# This is early Illumina data, I'm assuming qualities are solexa ASCII(QSolexa+64) scores. but they might not be. 
# they *might* be standard fastq 
task :run_bowtie, :roles => group_name do

  ip_i = ip.sub('.export.txt', '.fastq')
  input_i = input.sub('.export.txt', '.fastq')

  ip_o = ip.sub('.export.txt', '.sam')
  input_o = input.sub('.export.txt', '.sam')
  
  run "#{working_dir}/bowtie-0.12.7/bowtie --sam --al --best --solexa-quals  -q mm9 #{ip_i} > #{ip_o}"
  run "#{working_dir}/bowtie-0.12.7/bowtie --sam --al --best --solexa-quals  -q mm9 #{input_i} > #{input_o}"

end
before "run_bowtie", "EC2:start"



# fetch samtools from svn
desc "get samtools"
task :get_samtools, :roles => group_name do
  sudo "apt-get -y install subversion"
  run "svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
end
before "get_samtools", "EC2:start"


desc "build samtools"
task :build_samtools, :roles => group_name do
  sudo "apt-get -y install zlib1g-dev libncurses5-dev"
  run "cd /home/ubuntu/samtools && make"
end
before "build_samtools", "EC2:start"


desc "install samtools"
task :install_samtools, :roles => group_name do
  sudo "cp /home/ubuntu/samtools/samtools /usr/local/bin/samtools"
end
before "install_samtools", "EC2:start"


desc "make bam from sam"
task :to_bam, :roles => group_name do
  run "wget -O #{working_dir}/mm9_lengths  'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths'"
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.sam/)}
  files.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/mm9_lengths -o #{mount_point}/#{f_out} #{mount_point}/#{f}"
  }
end
before "to_bam", "EC2:start"



desc "sort bam"
task :sort_bam, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '_sorted')
    run "samtools sort #{mount_point}/#{f}  #{mount_point}/#{f_out}"
  }
end
before "sort_bam", "EC2:start"



desc "remove duplicates"
task :rmdups, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted\.bam/)}
  files.each{|f| 
    f_out = f.sub('sorted', 'sorted_nodups')
    run "cd #{mount_point} && samtools rmdup -s #{f}  #{f_out}"
  }
end
before "rmdups", "EC2:start"



desc "index bam files"
task :index, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '.bai')
    run "cd #{mount_point} && samtools index  #{f} #{f_out}"
  }
end
before "index", "EC2:start"



desc "download bam files"
task :get_bam, :roles => group_name do
  `rm -Rf results/alignment/bowtie` #remove previous results
  `mkdir -p results/alignment/bowtie`
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups/)}
  files.each{|f|
    download( "#{mount_point}/#{f}", "results/alignment/bowtie/#{f}")
  }
end
before "get_bam", 'EC2:start'


# Need to do this or macs will run out of space to make WIG files on a 30GB partition.
# If you want the intermediate files, use a bigger EBS vol
desc "cleanup intermediate files"
task :bam_tidy, :roles =>group_name do
  files = capture "ls #{mount_point}"
  
  #don't delete the macs results, original export.txt files or final bam,bai files
  files = files.split("\n")
  files = files.reject{
    |f| f.match(/^macs/) ||  
    f.match(/sorted_nodups/) || 
    f.match(/export.txt$/) }
  files.each{|f| run("rm -Rf #{f}")}
end 
before "bam_tidy", 'EC2:start'



### Macs ?

macs_url ="http://liulab.dfci.harvard.edu/MACS/src/MACS-1.4.0beta.tar.gz"
macs_version = "MACS-1.4.0beta"

task :install_macs, :roles => group_name do
  sudo "apt-get install -y python"
  run "cd #{working_dir} && wget --http-user macs --http-passwd chipseq #{macs_url}"
  run "cd #{working_dir} && tar -xvzf #{macs_version}.tar.gz"
  run "cd #{working_dir}/#{macs_version} && sudo python setup.py install"
  sudo "ln -s /usr/local/bin/macs* /usr/local/bin/macs"
end
before "install_macs", 'EC2:start'

task :install_peaksplitter, :roles => group_name do
  url ='http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp_1.0.tar.gz'
  filename = 'PeakSplitter_Cpp_1.0.tar.gz'
  bin = 'PeakSplitter_Cpp/PeakSplitter_Linux64/PeakSplitter'
  run "cd #{working_dir} && curl #{url} > #{filename}"
  run "cd #{working_dir} && tar -xvzf #{filename}"
  run "sudo cp #{working_dir}/#{bin} /usr/local/bin/PeakSplitter"
end 
before 'install_peaksplitter', 'EC2:start'

#you'll need to have done "install_r" and install_peak_splitter to do this
task :run_macs, :roles => group_name do

  treatment = "#{mount_point}/CMN039_sorted_nodups.bam"
  control = "#{mount_point}/CMN038_sorted_nodups.bam"
  genome = 'mm'
  bws = [300]
  pvalues = [0.00001]

  #unsure what p values and bandwidths are appropriate, try a few?
  bws.each {|bw|
    pvalues.each { |pvalue|

      dir = "#{mount_point}/macs_#{bw}_#{pvalue}"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --bw #{bw} --pvalue #{pvalue}"
      run "cd #{dir} && #{macs_cmd}"
      
      dir = "#{mount_point}/macs_#{bw}_#{pvalue}_subpeaks"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      # With SubPeak finding
      # this will take a lot longer as you have to save the wig file 
      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --call-subpeaks  --bw #{bw} --pvalue #{pvalue} --wig"
      run "cd #{dir} && #{macs_cmd}"

    }
  }
  
end
before 'run_macs', 'EC2:start'

#pack up the runs and downloads them to the server (without the wig files)
task :pack_macs, :roles => group_name do
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  macs_dirs.each{|d|
    run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
  }
  
end
before 'pack_macs','EC2:start' 

task :get_macs, :roles => group_name do
  macs_files = capture "ls #{mount_point}"
  macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  res_dir = 'results/alignment/bowtie/peakfinding/macs'
  `rm -Rf #{res_dir}`
  `mkdir -p #{res_dir}`
  macs_files.each{|f| 
    download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
    `cd #{res_dir} && tar -xvzf #{f}`
  }

end
before 'get_macs', 'EC2:start'



#if you want to keep the results

#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




