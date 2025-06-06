# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update c3g_tools module version
vim ../genpipes/resources/modules/c3g_tools.sh
(VERSION=2.12.8)

# Tag the branch and push the tag. 
git tag 2.12.8 -m 'Release 2.12.8'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 2.12.8"
git push

# Create a release tarball archive
git archive --format=tar --prefix=c3g_tools-2.12.8/ <latest_commmit> | gzip > ~/c3g_tools-2.12.8.tar.gz

# Create a release on GitHub, pointing to the new tag
# Include the archive in the release

# Version bump the value. Until the next release, add '-beta' e.g. 2.12.9-beta
vim VERSION
git commit -m "Version bump to 2.12.9-beta" VERSION
git push

# Deploy c3g_tools-<VERSION> as a module on all clusters

