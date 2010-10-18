require 'catpaws'
require 'find'
require 'fileutils'

require 'right_aws'

###################################################################
# Align the data to mm9
#
# CMN038 = NS5 Input
# CMN039 = NS5 H3K4me3 ChIP
#

#### config for catpaws ####

set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", }
set :nhosts, 1
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
#set :ami, 'emi-0243110F' #64bit Ubuntu 9.04 on Oxford NGS Eucalyptus. Any Ubuntu should be ok. 
set :ami, 'ami-52794c26' #32-bit ubuntu, EC2 eu-west-1
set :instance_type, 'm1.small'
set :group_name, 'chipseq'
set :working_dir, '/mnt/chipseq'
set :script_dir, working_dir+'/scripts'
role(:master, 'localhost')

set :s3_url, ENV['S3_URL']
set :s3_port, ENV['S3_PORT']
set :s3_protocol, ENV['S3_PROTOCOL']


s3 = RightAws::S3.new(aws_access_key, aws_secret_access_key,
                :server => s3_server,
		:port => s3_port,
		::protocol => s3_protocol
)

local_dir = ENV['DATADIR']+'/ngs_h3k4me3_chipseq'



## export.txt files are on s3, one file per bucket.
#buckets = AWS::S3::Bucket.list()
#bucket_names = buckets. map {|b| b.name}.grep(/^buckley\.gis100710.*/)
#object_names = bucket_names.map {|b| AWS::S3::Bucket.find(b).objects[0].key}
#
## we'll process each file on a sep instance
#set :nhosts, object_names.length
#
## s3 setup
#desc "install s3 client"
#task :install_s3, :roles => 'chipseq' do
#  sudo 'apt-get update'
#  sudo 'apt-get install -y s3cmd'
#end
#before 'install_s3', 'EC2:start'
#
#desc "upload the s3 config file"
#task :s3_config, :roles => 'chipseq' do
#  upload('/space/cassj/GIS_15_07_10/ChIPseq/s3cfg', '/home/ubuntu/.s3cfg')
#end 
#before 's3_config', 'EC2:start'
#
#
#
## fetch the data files
#task :fetch_files_setup, :roles => 'chipseq' do 
#  servers = find_servers_for_task(current_task)
#  it = 0..(servers.length - 1)
#  it.each do |i|
#    host = servers[i]
#    bucket = bucket_names[i]
#    object = object_names[i]
#    filename = bucket
#
#    cmd  = "s3cmd get s3://#{bucket}/#{object} #{working_dir}/#{object}"
#    put(cmd, "#{working_dir}/s3_fetch.sh", :hosts => host)
#    run("chmod +x #{working_dir}/s3_fetch.sh", :hosts => host)
#    
#    # to keep track of which starting file we're working from
#    run("export S3FILENAME=#{filename}")
#  end
#end
#before "fetch_files_setup", "EC2:start" 
#
#
#desc "Fetch an export file for each ec2 instance"
#task :fetch_files, :roles => 'chipseq' do 
#  run "#{working_dir}/s3_fetch.sh"
#end
#before "fetch_files", "fetch_files_setup"
#
#
## get bowtie
#task :get_bowtie, :roles => 'chipseq' do
#  run "wget 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.5/bowtie-0.12.5-linux-i386.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie%2F&ts=1282910294&use_mirror=kent' -O #{working_dir}/bowtie.zip"
#  run "cd #{working_dir} && unzip bowtie.zip"
#end 
#before 'get_bowtie', 'EC2:start'
#
#
## and the indexed genome
#
#task :fetch_genome, :roles => 'chipseq' do
#  run "s3cmd get --force s3://bowtie.mm9/mm9.ebwt.zip #{working_dir}/bowtie-0.12.5/indexes/mm9.ebwt.zip"
#  sudo "apt-get -y install unzip"
#  run "rm  #{working_dir}/bowtie-0.12.5/indexes/chr*"
#  run "cd  #{working_dir}/bowtie-0.12.5/indexes && unzip -o mm9.ebwt.zip"
#end 
#before "fetch_genome","EC2:start"
#
#
## convert export.txt files to fastq for bowtie.
#task :export2fastq, :roles => 'chipseq' do
#  run "wget 'http://github.com/cassj/my_bioinfo_scripts/raw/master/ngs/illumina/export2fastq.pl' -O #{working_dir}/export2fastq.pl"
#  run "chmod +x #{working_dir}/export2fastq.pl"
#  run "#{working_dir}/export2fastq.pl #{working_dir}/*export.txt > #{working_dir}/export.fastq"
#end 
#before "export2fastq", "EC2:start"
#
#
#
## run bowtie on the fastq file
#task :do_bowtie, :roles => 'chipseq' do
#  run "#{working_dir}/bowtie-0.12.5/bowtie --sam --al --best --phred64-quals -q mm9 #{working_dir}/export.fastq > #{working_dir}/export.sam"
#end 
#before "do_bowtie", "EC2:start"
#
#
##### I am HERE! ####
#
#
#
#
## fetch samtools from svn
#desc "get samtools"
#task :get_samtools, :roles => :chipseq do
#  sudo "apt-get -y install subversion"
#  run "svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
#end
#before "get_samtools", "EC2:start"
#
#
#desc "build samtools"
#task :build_samtools, :roles => :chipseq do
#  sudo "apt-get -y install zlib1g-dev libncurses5-dev"
#  run "cd /home/ubuntu/samtools && make"
#end 
#before "build_samtools", "EC2:start"
#
#
#desc "install samtools"
#task :install_samtools, :roles => :chipseq do
#  sudo "cp /home/ubuntu/samtools/samtools /usr/local/bin/samtools"
#end
#before "install_samtools", "EC2:start"
#
#
#desc "make bam from sam"
#task :to_bam, :roles => :chipseq do
#  run "wget -O #{working_dir}/mm9_lengths  'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths'"
#  run "samtools view -bt #{working_dir}/mm9_lengths -o #{working_dir}/export.bam #{working_dir}/export.sam"
#end 
#before "to_bam", "EC2:start"
#
#desc "sort bam"
#task :sort_bam, :roles => :chipseq do
#  run "samtools sort #{working_dir}/export.bam #{working_dir}/export_sorted"
#end 
#before "sort_bam", "EC2:start"
#
#desc "remove duplicates"
#task :rmdups, :roles => :chipseq do
#  run "samtools rmdup -s #{working_dir}/export_sorted.bam #{working_dir}/export_nodups.bam"
#end 
#before "rmdups", "EC2:start"
#
#desc "index bam files"
#task :index, :roles => :chipseq do
#  run "samtools index #{working_dir}/export_nodups.bam #{working_dir}/export_nodups.bai" 
#end 
#before "index", "EC2:start"
#
#
#desc "upload to s3"
#task "to_s3", :roles => :chipseq do
#
#  servers = find_servers_for_task(current_task)
#  puts servers
#  it = 0..(servers.length - 1)
#  it.each do |i|
#    host = servers[i]
#    bucket = bucket_names[i]
#    object = object_names[i]
#    bucket = "bam."+bucket
#    run("s3cmd mb s3://#{bucket}", :hosts => host)
#    run("s3cmd put #{working_dir}/export.bam s3://#{bucket}", :hosts => host)
#    run("s3cmd put #{working_dir}/export_sorted.bam s3://#{bucket}", :hosts => host)
#    run("s3cmd put #{working_dir}/export_nodups.bam s3://#{bucket}", :hosts => host)
#    run("s3cmd put #{working_dir}/export_nodups.bai s3://#{bucket}", :hosts => host)
#  end
#
#end
#before "to_s3", "EC2:start"
#
#
#
#
#
#


