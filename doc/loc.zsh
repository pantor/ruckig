echo "Lines of Code of Ruckig"
cloc .. --exclude-lang=HTML,CSS,JSON,XML,YAML,zsh,Markdown --exclude-dir=reflexxes,build,examples,test

echo
echo "Lines of Code of Reflexxes"
cloc ../reflexxes --exclude-lang=HTML,CSS,JSON,XML,YAML,zsh,Markdown --fullpath --not-match-d=ReflexxesTypeIV/Windows_VS2010Express,ReflexxesTypeIV/src/RMLPositionSampleApplications,ReflexxesTypeIV/src/RMLVelocitySampleApplications,ReflexxesTypeIV/src/TypeIVRMLVelocity
