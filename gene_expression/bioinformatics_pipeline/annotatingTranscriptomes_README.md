# code I ran to start the annotating transcriptomes code

cd programs/
wget https://github.com/z0on/annotatingTranscriptomes/archive/master.zip
unzip master.zip

# pwd = /scratch/projects/and_transcriptomics/programs/annotatingTranscriptomes-master

# first pl script to run
seq_stats.pl
# first few lines show perl modules i don't have

#path to fasta file:
/scratch/projects/and_transcriptomics/programs/annotatingTranscriptomes-master/seq_stats.pl Acropora_cervicornis.mrna-transcripts.fa

# need to install a couple perl modules (i.e. Bio::Perl) but i don't think I can do that with the shared version of perl from Pegasus,
#so I'm going to try installing perl locally in my scratch space

wget https://www.cpan.org/src/5.0/perl-5.32.1.tar.gz
tar -xzvf perl-5.32.1.tar.gz
cd perl-5.32.1
./Configure -des -Dprefix=$HOME/localperl
make
make test
make install
nano ~/.bash_profile
export PATH="$HOME/localperl/bin:$PATH"
source ~/.bash_profile
perl -v
which perl # this should now be the locally installed one

# now install cpanminus for installing Bio::perl

cpan App::cpanminus
cpanm --local-lib=~/perl5 local::lib
eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
echo 'eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)' >> ~/.bash_profile
source ~/.bash_profile

cpanm Bio::Perl

# ok none of this worked. I think I can't run this on Pegasus.
