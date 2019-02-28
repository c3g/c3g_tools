# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_tools module version
vim ../genpipes/resources/modules/mugqic_tools.sh
(VERSION=2.2.2)

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 2.2.2 -m 'Release 2.2.2'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 2.2.2"

# Create a release tarball archive
git archive --format=tar --prefix=mugqic_tools-2.2.2/ 2.2.2 | gzip > ~/mugqic_tools-2.2.2.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/mugqic_tools/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 2.2.3-beta
vim VERSION
git commit -m "Version bump to 2.2.3-beta" VERSION
git push

# Deploy mugqic_tools-<VERSION> as a module on all clusters

# Send a message to the mailing list:
mugqic_pipelines@googlegroups.com

# In JIRA, add a release date to the 'Version' category of the administer project
