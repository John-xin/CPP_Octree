# CPP_Octree

## download git
1. https://git-scm.com/download/win
2. select "64-bit Git for Windows Portable"
3. install

## create ssh key
1. input this command and follow instructions\
`ssh-keygen -t rsa -C "ewangxin11@gmail.com"`

2. add the key to github account\
https://github.com/settings/keys

3. check if successfully added by this command\
`ssh -T git@github.com`

## create a repository in github
1. new a repository
2. copy the ssh address not http

## use git to upload code
1. set global user name and email\
```
git config --global user.name 这里是你的用户名
git config --global user.email example@email.com
```
2. enter code folder and init a git repository\
```
cd E:
cd codeFolder
git init //初始化为git仓库
```
.git should be created under codeFolder
  
3a. clone code from github repository to current folder\
`git clone -b <branchname> <remote-repo-url> // replace this address`


or

3b. upload my code to github repository 
```
git add . //压入缓存，或者说工作区
git commit -m "first commit" //提交入管理器，或者说暂存区 
git remote add origin aSshAddress //create name representing this long address
git push -u origin master //以主干master提交代码>  
``` 
  
## useful git commands
```
git remote -v //list all names
git remote remove aName //delete name
```
 
