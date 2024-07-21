# Group3_Final_Report

This is a analysis of RNA-seq, metabolic with respect to metastasis of cancer.

最後的分析報告放在 \~/Analysis/Final Report.html

## 資料夾敘述

這個專案中有proposal和analysis兩個資料夾。兩個資料夾中Proposal是我們針對資料分析所提出來的proposal，實作的分析都放在analysis中。關於分析的qmd檔如何執行會在下一節中敘述。

## 資料分析

在analysis的資料中是由一份主要的qmd檔去讀取每個組員各自做的qmd檔。主要的那份qmd檔案是叫做 *final report.qmd* ，裡面讀取了六份qmd檔案，不用執行各自的檔案，僅需執行*final report.qmd* 即可完成整份檔案的編譯。

此外，此份文件使用了r和python的兩種語言，因此需要指定python的虛擬環境才可執行檔案，詳情在下一節中敘述。

## 關於python的執行

要同時執行python和r需要使用到r語言中`reticulate`這個套件，因此在*final report.qmd*的檔案中，一開始有呼叫這個套件。這個套件可以指定一個電腦中存在的python虛擬環境，並在環境中安中所需要python套件即可（虛擬環境的位置需要更改以執行檔案）。
