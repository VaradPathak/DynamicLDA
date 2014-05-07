DynamicLDA
==========

Dynamic Topic Model of Reuters News Articles between 2007-2013
--------------------------------------------------------------
<p>We have implemented fast version of <a href="https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDAQFjAA&url=http%3A%2F%2Fwww.cs.cmu.edu%2F~lafferty%2Fpub%2Fdtm.pdf&ei=YZJqU_-ABIL_oQTVroCQDg&usg=AFQjCNGicS7Nr_Q76R5uSUczaUP2DaAd1A&sig2=xOoJWejgXXVBTu9wf4vAVw&bvm=bv.66111022,d.cGU&cad=rja">Dynamic Topic Model</a> proposed by David Blei and John Lafferty in 2006.</p>
<p>This version takes advantage of new advancements in LDA model. We have implemented the LDA part using <a href="http://arxiv.org/abs/1305.2452">SCVB0</a> proposed by Foulds, et al 2013. This is parallelized implementation of SCVB0 using OpenMP.</p>
<p>As per our evaluation, even our Serial version gives 36X speedup and the Parallel version when run on core 2 duo machine 2GHz 2Gb gives 53X speedup.</p>
<p> (<a href="https://www.dropbox.com/s/8hudevuubex7egq/DTM%20Final%20Report.pdf">Report with detail evaluation</a>)</p>

Reuters News Dataset Details
----------------------------
Timestamped News articles published by Reuters between 2007 and 2013. This is corpus of 161,989 documents with vocab size of 32,468 after preprocessing. Following are the preprocessing steps performed (Scripts are available in Scrapper folder)

 - From Reuters data we removed all the docs which have length less than 100 words
 - We have scrapped random 10% of the data from each day. This was done just to minimize the corpus size.The assumption is that randomly selected data wont cause problem while find the long and big topics.
 - We removed all the punctuation marks and performed stemming using Porter2 stemmer
 - We also removed the words which have frequency of less tan 25 or more than 100,000
example run of text2ldac:

Topic Chains
------------
We have investigated the <a href="http://dl.acm.org/citation.cfm?id=1964765">Topic Chains</a> a solution to topic Birth-Death problem in Dynamic LDA proposed by Kim, at al in 2013
 - We use the same Reuters dataset and use the Jensen-Shannon (JS) divergence to compare similarity between the topics.
 - We investigate the performance at different similarity thresholds and Window sizes and find similar results as given in the original paper
 - We identify some issues in the method and propose solutions to the same (Please refer the report for more details)

Execution Commands
------------------
 - Scrape Data from reuters archive website between startMonth for num_of_months<br>
<code>python __init__.py startMonth num_of_months</code>
 - Get Stopwords <code>python removeInfrequentWords.py</code>
 - Convert the text data to ldac format used by <a href="https://code.google.com/p/princeton-statistical-learning/downloads/detail?name=dtm_release-0.8.tgz">Blei's implementation</a><br>
<code>python multitext2ldac.py data_folder --stopwords stopwords_file</code>
 - Convert data to <a href="https://archive.ics.uci.edu/ml/datasets/Bag+of+Words">UCI</a> format
  <code>python ldac2uci.py</code>
 - Compile Dynamic LDA. <code>make</code>
 - Execute Dynamic Topic Modeling on UCI dataset<br>
 <code>./fastLDA UCIFormat_data_file iterations NumOfTopics MiniBatchSize Vocab_file GeneratePi</code>
 - Get the word trend in a topic<br>
 <code>python getWordVariation.py TopicId WordId PiFolderPath StartYear EndYear</code>
 - Compile Topic Chains GetData to get all the Topics in the dataset for all the TimeSlices
 <code>make GetData</code>
 - Execute GetData for Topic Chains<br> <code>./GetData UCIFormat_data_file iterations NumOfTopics MiniBatchSize Vocab_file GeneratePi</code>
 - Compile GenerateChains for Topic Chains <code>make GenerateChains</code>
 - Execute GenerateChains <code>./GenerateChains  Pi_folder num_topics WindowSize SimilarityThreshold</code>
