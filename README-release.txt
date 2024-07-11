# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_tools module version
vim ../genpipes/resources/modules/mugqic_tools.sh
(VERSION=2.12.8)

# Tag the branch and push the tag. 
git tag 2.12.8 -m 'Release 2.12.8'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 2.12.8"
git push

# Create a release tarball archive
git archive --format=tar --prefix=mugqic_tools-2.12.8/ <latest_commmit> | gzip > ~/mugqic_tools-2.12.8.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/mugqic_tools/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 2.12.9-beta
vim VERSION
git commit -m "Version bump to 2.12.9-beta" VERSION
git push

# Deploy mugqic_tools-<VERSION> as a module on all clusters

