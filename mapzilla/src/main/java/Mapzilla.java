package main.java;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.Socket;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Mapzilla {
	
	public static final int BLOCK_SIZE = 2000;
	public static final int READ_LENGTH = 100;
	public static final int CHUNK_SIZE = 25;
	public static final long READS_PENDING_OUTPUT_MAX = 1000000;
	
	public static String HOSTNAME = "hm1-pvt.bioinf.unc.edu";
//	public static String HOSTNAME = "localhost";
	
//	public static final int SERVER_PORT = 5000;
	public static final int SERVER_PORT = 6000;
	
	private OutputRunnable outputThread;
		
//	int numThreads = 10;
	
//	int numThreads = 1;
	
	int numThreads = 14;
	
//	private Map<String, byte[]> refMap = new HashMap<String, byte[]>();
//	private String refFileName = "/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/hg19.fa";
	
	private List<RefPosition> refPositions = new ArrayList<RefPosition>();
	
	private long startTime;
	
	private void initPositions() {
		initPos("chr1", 249250621);
		initPos("chr2", 243199373);
		initPos("chr3", 198022430);
		initPos("chr4", 191154276);
		initPos("chr5", 180915260);
		initPos("chr6", 171115067);
		initPos("chr7", 159138663);
		initPos("chr8", 146364022);
		initPos("chr9", 141213431);
		initPos("chr10", 135534747);
		initPos("chr11", 135006516);
		initPos("chr12", 133851895);
		initPos("chr13", 115169878);
		initPos("chr14", 107349540);
		initPos("chr15", 102531392);
		initPos("chr16", 90354753);
		initPos("chr17", 81195210);
		initPos("chr18", 78077248);
		initPos("chr19", 59128983);
		initPos("chr20", 63025520);
		initPos("chr21", 48129895);
		initPos("chr22", 51304566);
		initPos("chrX", 155270560);
		initPos("chrY", 59373566);
		initPos("chrM", 16571);
	}
	
	/*
	private void initPositions() {
		initPos("chr7", 159138663);
	}
	*/
	
	private void initPos(String ref, int length) {
		long offset = 0;
		for (RefPosition refPos : refPositions) {
			offset += refPos.length;
		}
		
		RefPosition refPos = new RefPosition();
		refPos.ref = ref;
		refPos.length = length;
		refPos.offset = offset;
		
		refPositions.add(refPos);
	}
	
	private RefPosition findRefPosition(long pos) {
		for (RefPosition refPos : refPositions) {
			if ((pos >= refPos.offset) && (pos < refPos.offset + refPos.length)) {
				return refPos;
			}
		}

		return null;
		//throw new IllegalArgumentException("Invalid reference position: " + pos);
	}
	
	private void initSockets(Socket[] sockets, OutputStream[] os, BufferedOutputStream[] bos,
			DataOutputStream[] dos, InputStream[] is, BufferedInputStream[] bis, DataInputStream[] dis) throws IOException {
		
		for (int i=0; i<numThreads; i++) {
			sockets[i] = new Socket(HOSTNAME, SERVER_PORT);
			os[i] = sockets[i].getOutputStream();
			bos[i] = new BufferedOutputStream(os[i]);
			dos[i] = new DataOutputStream(bos[i]);
			is[i] = sockets[i].getInputStream();
			bis[i] = new BufferedInputStream(is[i]);
			dis[i] = new DataInputStream(bis[i]);
		}
	}
	
	private void cleanupSockets(Socket[] sockets, OutputStream[] os, BufferedOutputStream[] bos,
			DataOutputStream[] dos, InputStream[] is, BufferedInputStream[] bis, DataInputStream[] dis) throws IOException {
		
		for (int i=0; i<numThreads; i++) {
			
			dos[i].writeInt(-1);
			
			dos[i].close();
			bos[i].close();
			os[i].close();
			dis[i].close();
			bis[i].close();
			is[i].close();
			sockets[i].close();
		}
	}
	
	private long elapsedSecs() {
		return (System.currentTimeMillis() - startTime) / 1000;
	}
	
	public void run(String input) throws IOException, InterruptedException {
		
		startTime = System.currentTimeMillis();
		
		initPositions();
		
		Socket[] sockets = new Socket[numThreads];
		OutputStream[] os = new OutputStream[numThreads];
		BufferedOutputStream[] bos = new BufferedOutputStream[numThreads];
		DataOutputStream[] dos = new DataOutputStream[numThreads];
		
		InputStream[] is = new InputStream[numThreads];
		BufferedInputStream[] bis = new BufferedInputStream[numThreads];
		DataInputStream[] dis = new DataInputStream[numThreads];
		
		initSockets(sockets, os, bos, dos, is, bis, dis);
		
		outputHeader();
		// Open input file
//		BufferedReader reader = new BufferedReader(new FileReader(input));
		FastqInputFile fastq = new FastqInputFile();
		fastq.init(input);
		
		outputThread = new OutputRunnable();
		Thread thread = new Thread(outputThread);
		thread.start();
		
		LookupRunnable[] lookups = new LookupRunnable[numThreads];
		Thread[] lookupThreads = new Thread[numThreads];
		
		for (int i=0; i<numThreads; i++) {
			lookups[i] = new LookupRunnable(dis[i], dos[i], this);
			lookupThreads[i] = new Thread(lookups[i]);
			lookupThreads[i].start();
		}
				
		List<FastqRecord> reads = new ArrayList<FastqRecord>();
		FastqRecord read = fastq.getNextRecord();
//		String line = reader.readLine(); 
		int i = 0;
		
		while (read != null) {
			reads.add(read);
			
			if (reads.size() == BLOCK_SIZE) {
//				lookup(reads, dos, dis);
				
				int idx = 0;
				boolean isPending = true;
				boolean isWaitingOnOutput = false;
				
				while (isPending) {
					
					if (outputThread.getPendingCount() >= READS_PENDING_OUTPUT_MAX) {
						if (!isWaitingOnOutput) {
							isWaitingOnOutput = true;
							System.err.println("Waiting for output thread to catch up: " + elapsedSecs());
						}
						sleep(100);
					} else {
						
						if (isWaitingOnOutput) {
							isWaitingOnOutput = false;
							System.err.println("Done waiting for output thread to catch up: " + elapsedSecs());
						}
						
						if (!lookups[idx].hasNewReadsQueued()) {
							lookups[idx].setNewReads(reads);
							isPending = false;
						} else {
							if (idx == numThreads-1) {
								sleep(10);
								idx = 0;
							} else {
								idx++;
							}
						}
					}
				}
								
//				lookup1.setNewReads(reads);
				i++;
				reads = new ArrayList<FastqRecord>();  // Create new reference to avoid stomping on worker threads
			}
//			line = reader.readLine();
			read = fastq.getNextRecord();
		}
		
		// Send leftovers
		int idx = i % numThreads;
		lookups[idx].setNewReads(reads);

		for (int j=0; j<numThreads; j++) {
			lookups[j].done();
		}
		
		System.err.println("Processing done...");
		
//		reader.close();
		fastq.close();
		
		for (int j=0; j<numThreads; j++) {
			System.err.println("Waiting for lookup thread: " + j);
			lookupThreads[j].join();
		}
		
		outputThread.done();
		
		System.err.println("Sending done notification to server");
		// Indicate that we are done
		cleanupSockets(sockets, os, bos, dos, is, bis, dis);
		
		System.err.println("Waiting for output thread");
		thread.join();
		System.err.println("Done");
	}
	
	private void outputHeader() {
		System.out.println("@SQ	SN:chr1	LN:249250621");
		System.out.println("@SQ	SN:chr2	LN:243199373");
		System.out.println("@SQ	SN:chr3	LN:198022430");
		System.out.println("@SQ	SN:chr4	LN:191154276");
		System.out.println("@SQ	SN:chr5	LN:180915260");
		System.out.println("@SQ	SN:chr6	LN:171115067");
		System.out.println("@SQ	SN:chr7	LN:159138663");
		System.out.println("@SQ	SN:chr8	LN:146364022");
		System.out.println("@SQ	SN:chr9	LN:141213431");
		System.out.println("@SQ	SN:chr10	LN:135534747");
		System.out.println("@SQ	SN:chr11	LN:135006516");
		System.out.println("@SQ	SN:chr12	LN:133851895");
		System.out.println("@SQ	SN:chr13	LN:115169878");
		System.out.println("@SQ	SN:chr14	LN:107349540");
		System.out.println("@SQ	SN:chr15	LN:102531392");
		System.out.println("@SQ	SN:chr16	LN:90354753");
		System.out.println("@SQ	SN:chr17	LN:81195210");
		System.out.println("@SQ	SN:chr18	LN:78077248");
		System.out.println("@SQ	SN:chr19	LN:59128983");
		System.out.println("@SQ	SN:chr20	LN:63025520");
		System.out.println("@SQ	SN:chr21	LN:48129895");
		System.out.println("@SQ	SN:chr22	LN:51304566");
		System.out.println("@SQ	SN:chrX	LN:155270560");
		System.out.println("@SQ	SN:chrY	LN:59373566");
		System.out.println("@SQ	SN:chrM	LN:16571");
	}
	
	public long hash(byte[] bytes) {
		long hash = 0;
		
		for (byte b : bytes) {

			hash = hash << 2;
			long val = to2bit(b);
			hash += val;
		}
		
		return hash;
	}
	
	private long to2bit(byte b) {
		long val;
		
		switch (b) {
			case 'A':
				val = 0;
				break;
			case 'C':
				val = 1;
				break;
			case 'T':
				val = 2;
				break;
			case 'G':
				val = 3;
				break;
			default:
				throw new IllegalArgumentException("Invalid base: " + b);
		}

		return val;
	}

	
	private long[] getHashes(String sequence) {
		
		long[] hashes = new long[sequence.length() / CHUNK_SIZE];
		
		int idx = 0;
		for (int i=0; i<sequence.length(); i+=CHUNK_SIZE) {
			String str = sequence.substring(i, i+CHUNK_SIZE);
			
			if (str.contains("N")) {
				//TODO: Use more "special" number here?
				hashes[idx++] = 0;
			} else {			
				byte[] bytes = str.getBytes();
				hashes[idx++] = hash(bytes);
			}
		}
		
		return hashes;
	}
	
	private void lookup(List<FastqRecord> reads, DataOutputStream dos, DataInputStream dis) throws IOException {
//		System.err.println("Sending: " + reads.size()*4);
		dos.writeInt(reads.size() * 4);
		for (FastqRecord read : reads) {
			
			long[] hashes = getHashes(read.getSequence());
			for (int i=0; i<hashes.length; i++) {
				dos.writeLong(hashes[i]);
			}
			
//			dos.writeBytes(read.getSequence());
		}
		dos.flush();
		
//		System.err.println("Sent...");
		
//		int[] positions = new int[reads.size() * 4];
		
		long[] positions = new long[4];
		int idx = 0;
		int numPositions = reads.size() * 4;
		
		List<String> outputReads = new ArrayList<String>();
		
		for (int i=0; i<=numPositions; i++) {

			if (idx == 4) {
				long position = -1;
				for (idx=0; idx<positions.length-1 && position==-1; idx++) {
					
					for (int idx2=idx+1; idx2<positions.length; idx2++) {
						if (positions[idx] == positions[idx2] - (idx2-idx)*CHUNK_SIZE) {
							position = positions[idx] - (idx * CHUNK_SIZE) + 1;
							break;
						}
					}
					
//					System.out.println("position: " + positions[idx]);
				}
				FastqRecord read = reads.get((i-1)/4);
				
				String readToOutput = getReadToOutput(read, position);
				outputReads.add(readToOutput);
				
//				if (outputReads.size() == BLOCK_SIZE) {
//					outputThread.setNewReadsToOutput(outputReads);
//					outputReads = new ArrayList<String>();
//				}
				idx = 0;
			}
			
			if (i < numPositions) {
				long s = System.currentTimeMillis();
				try {
					positions[idx++] = dis.readLong();
//					System.out.println("pos: " + positions[idx-1]);
				} catch (EOFException x) {
					long e = System.currentTimeMillis();
					System.out.println("wait time: " + (e-s));
					throw x;
				}
//				System.err.println("recieved: " + positions[idx - 1]);
			}			
		}
		
//		System.err.println("Finished reading response");
		
		if (outputReads.size() > 0) {
			outputThread.setNewReadsToOutput(outputReads);
		}
	}
	
	/*
	private void lookup(List<FastqRecord> reads, DataOutputStream dos, DataInputStream dis) throws IOException {
		dos.writeInt(reads.size());
		for (FastqRecord read : reads) {
			dos.writeBytes(read.getSequence());
		}
		dos.flush();
		
//		int[] positions = new int[reads.size() * 4];
		
		// # reads * 4 chunks per read * 4 bytes per chunk
		int totalBytes = reads.size() * 4 * 4;
		byte[] bytes = new byte[totalBytes];
		
		int offset = 0;
		int bytesRead = dis.read(bytes, offset, totalBytes);
		if (bytesRead == -1) {
			throw new RuntimeException("Premature EOF.  bytesRead: " + bytesRead);
		}
		while (bytesRead < totalBytes) {
			offset = bytesRead;
			int readBytes = dis.read(bytes, offset, totalBytes - bytesRead);
			if (readBytes == -1) {
				throw new RuntimeException("Premature EOF.  bytesRead: " + bytesRead);
			}
			bytesRead += readBytes;
		}

		ByteBuffer byteBuf = ByteBuffer.wrap(bytes);
		IntBuffer intBuf = byteBuf.asIntBuffer();
		int[] positions = new int[intBuf.remaining()];
		intBuf.get(positions);
		
//		for (int i=0; i<positions.length; i++) {
//			System.out.println("pos: " + i + " - " + positions[i]);
//		}
		
		List<String> outputReads = new ArrayList<String>();
		
		int chunkIdx = 0;
		for (int i=0; i<=positions.length; i++) {

			if (chunkIdx == 4) {
				int position = -1;
				
				int chunkStart = i-4;
				int chunkStop  = i;
				
				for (int idx=chunkStart; idx<chunkStop; idx++) {
					
					for (int idx2=idx+1; idx2<chunkStop; idx2++) {
						if (positions[idx] == positions[idx2] - (idx2-idx)*CHUNK_SIZE) {
							position = positions[idx] - (idx * CHUNK_SIZE) + 1;
							break;
						}
					}
				}
				FastqRecord read = reads.get((i-1)/4);
				String readToOutput = getReadToOutput(read, position);
				outputReads.add(readToOutput);
				if (outputReads.size() == BLOCK_SIZE) {
					outputThread.setNewReadsToOutput(outputReads);
					outputReads = new ArrayList<String>();
				}
				
				chunkIdx = 0;
			}
			
			chunkIdx++;
		}
		
		if (outputReads.size() > 0) {
			outputThread.setNewReadsToOutput(outputReads);
		}
	}
	*/
	
	private String getReadToOutput(FastqRecord read, long position) {
		
		RefPosition refPos = null;
		long pos = 0;
		if (position > 0) {
			refPos = findRefPosition(position);
			if (refPos != null) {
				pos = position - refPos.offset;
			} else {
				position = -1;
			}
		}

//		System.out.println("READ\t" + read.getId() + "\t7\t" + position + "\t" + read.getSequence());
		StringBuffer str = new StringBuffer(300);
		str.append(read.getId().substring(1, read.getId().length()));
		str.append('\t');
		str.append(position > 0 ? 0 : 4);
		str.append('\t');
		//str.append("chr7");
		str.append(refPos != null ? refPos.ref : "*");
		str.append('\t');
		str.append(pos > 0 ? pos : "*");
		str.append("\t50\t100M\t*\t*\t*\t");
		str.append(read.getSequence());
		str.append('\t');
		str.append(read.getQuality());
		return str.toString();
	}
	
	private static void sleep(int millis) {
		try {
			Thread.sleep(millis);
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	/*
	static class LookupRunnable implements Runnable {
		
		public void run() {
			
			while (!isDone()) {
		}
	}
	*/
	
	static class LookupRunnable implements Runnable {

		private List<FastqRecord> reads = null;
		private List<List<FastqRecord>> newReads = new ArrayList<List<FastqRecord>>();
		private boolean isDone = false;
		private DataInputStream dis;
		private DataOutputStream dos;
		private Mapzilla mapzilla;
		private static int i = 1;
		private int j = i++;
		
		public LookupRunnable(DataInputStream dis, DataOutputStream dos, Mapzilla mapzilla) {
			this.dis = dis;
			this.dos = dos;
			this.mapzilla = mapzilla;
//			System.out.println("j: " + j);
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			
			while (!isDone()) {
				reads = getNewReads();
				if (reads != null) {
					try {
//						System.out.println("j: " + j);
						mapzilla.lookup(reads, dos, dis);
					} catch (IOException e) {
						e.printStackTrace();
						throw new RuntimeException(e);
					}
					
					reads = null;
				} else {
					sleep(100);
				}
			}
		}
		
		public synchronized void setNewReads(List<FastqRecord> reads) {
			newReads.add(reads);
		}
		
		public synchronized List<FastqRecord> getNewReads() {
			if (newReads.size() > 0) {
				reads = newReads.remove(0);
			} else {
				reads = null;
			}
			
			return reads;
		}
		
		public synchronized boolean hasNewReadsQueued() {
			return newReads.size() > 0;
		}
		
		public synchronized boolean isDone() {
//			boolean _isDone = (newReadsToOutput == null) && isDone;
//			System.out.println("newReadsToOutput: " + newReadsToOutput);
//			System.out.println("isDone: " + isDone);
//			System.out.println("Done? " + _isDone);
//			return _isDone;
			
			return (newReads.size() == 0) && isDone;
		}
		
		public synchronized void done() {
			isDone = true;
		}
	}

	
	static class OutputRunnable implements Runnable {

		private List<String> readsToOutput = null;
		private List<String> newReadsToOutput = null;
		private boolean isDone = false;
		
		@Override
		public void run() {
			try {
//				BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
				
				int count = 0;
				while (!isDone()) {
					readsToOutput = getNewReadsToOutput();
					if (readsToOutput != null) {
						for (String read : readsToOutput) {
							System.out.println(read);
							if ((count % 250000) == 0) {
								System.out.flush();
							}
//							writer.write(read);
//							writer.write('\n');
							count++;
						}
						
						readsToOutput = null;
					} else {
						sleep(100);
					}
				}
				
//				writer.flush();
//				writer.close();
				
				if (1 == 2) {
					throw new IOException();
				}
				
				System.err.println("Reads output: " + count);
				System.err.println("New count: " + newCount);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		int newCount = 0;
		
		public synchronized void setNewReadsToOutput(List<String> reads) {
			newCount += reads.size();
			if (newReadsToOutput == null) {
				newReadsToOutput = reads;
			} else {
//				System.err.println("Adding reads!");
				newReadsToOutput.addAll(reads);
			}
		}
		
		public synchronized List<String> getNewReadsToOutput() {
			if (newReadsToOutput != null) {
				readsToOutput = newReadsToOutput;
				newReadsToOutput = null;
			}
			
			return readsToOutput;
		}
		
		public synchronized long getPendingCount() {
			return newReadsToOutput != null ? newReadsToOutput.size() : 0L;
		}
		
		public synchronized boolean isDone() {
//			boolean _isDone = (newReadsToOutput == null) && isDone;
//			System.out.println("newReadsToOutput: " + newReadsToOutput);
//			System.out.println("isDone: " + isDone);
//			System.out.println("Done? " + _isDone);
//			return _isDone;
			
			return (newReadsToOutput == null) && isDone;
		}
		
		public synchronized void done() {
			isDone = true;
		}
	}
	
	static class RefPosition {
		String ref;
		long offset;
		long length;
	}
	
	static class FastqInputFile {
	    
	    private static final int CACHING_DISABLED = -1;
	    
	    private BufferedReader reader;
	    private Map<Integer, FastqRecord> records = new HashMap<Integer, FastqRecord>();
	    private int recordNum = 0;
	    private int maxCachedLines;
	    
	    public void init(String filename, int maxCachedLines) throws FileNotFoundException {
	        openFile(filename);
	        this.maxCachedLines = maxCachedLines;
	    }
	    
	    public void init (String filename) throws FileNotFoundException {
	        init(filename, CACHING_DISABLED);
	    }
	    
	    private void openFile(String filename) throws FileNotFoundException {
	        reader = new BufferedReader(new FileReader(filename));
	    }

	    public FastqRecord getNextRecord() throws IOException {
	        String[] lines = new String[FastqRecord.NUM_LINES];
	        
	        for (int i=0; i<4; i++) {
	            String line = reader.readLine();
	            if (line == null) {
	                return null;
	            }
	            lines[i] = line;
	        }
	        
	        FastqRecord fastqRecord = new FastqRecord(lines);
	        return fastqRecord;
	    }
	    
	    /**
	     *  Retrieves the record 
	     */
	    public FastqRecord getRecord(int idx) throws IOException {
	        
	        if (isCachingDisabled()) {
	            throw new UnsupportedOperationException("Retrieving record by index not supported when caching is disabled.");
	        }
	        
	        if (idx < 1) {
	            throw new IllegalArgumentException("Read index cannot be less than 1.  Index: [" + idx + "]");
	        }
	        
	        FastqRecord record = null;
	        
	        // Attempt to get record from cache
	        if (records.containsKey(idx)) {
	            return records.get(idx);
	        }
	        
	        // Can't "read back" before the cache
	        if (idx <= recordNum) {
	            throw new UnsupportedOperationException("Cannot read back to an uncached record. Index: [" + idx + "]");
	        }
	        
	        // Read forward to the desired record and update the cache
	        while (idx > recordNum) {
	            record = getNextRecord();
	            recordNum++;
	            // Cache the newly read record and purge stale records.
	            // i.e. Put record 101 and remove record 1.  
	            records.put(recordNum, record);
	            records.remove(recordNum - maxCachedLines);
	        }
	        
	        // The cache should never be bigger than maxCachedLines
	        assert(records.size() <= maxCachedLines);
	        
	        return record;
	    }
	    
	    public void close() throws IOException {
	        if (reader != null) {
	            reader.close();
	        }
	    }
	    
	    private boolean isCachingDisabled() {
	        return maxCachedLines == CACHING_DISABLED;
	    }
	}

	static class FastqRecord {
	    
	    public static final int NUM_LINES = 4;
	    
	    private String[] lines = new String[NUM_LINES];
	    
	    public FastqRecord(String[] lines) {
	        if (lines.length != NUM_LINES) {
	            throw new IllegalArgumentException("Invalid number of lines for FastqRecord: [" + Arrays.toString(lines) + "]");
	        }
	        this.lines = lines;
	    }
	    
	    FastqRecord(String id, String sequence, String quality) {
	    	lines[0] = id;
	    	lines[1] = sequence;
	    	lines[2] = "+";
	    	lines[3] = quality;
	    }
	    
	    public String getId() {
	        return lines[0];
	    }
	    
	    public String[] getLines() {
	        return lines;
	    }

	    
	    /**
	     * Returns the portion of the id string leading up to "/"
	     */
	    public String getBaseId() {
	        int slashIdx = getId().indexOf("/");
	        int spaceIdx = getId().indexOf(" ");
	        
	        if ((slashIdx == -1) && (spaceIdx == -1)) {
	            return getId();
	        }
	        
	        int idx = -1;
	        if (slashIdx == -1) {
	            idx = spaceIdx;
	        } else if (spaceIdx == -1) {
	            idx = slashIdx;
	        } else {
	            idx = spaceIdx < slashIdx ? spaceIdx : slashIdx;
	        }
	        
	        return getId().substring(0, idx);
	    }
	    
	    /**
	     * Returns true if this FastqRecord has the same base id as the input FastqRecord 
	     */
	    public boolean hasSameBaseId(FastqRecord rec) {
	        return rec != null && this.getBaseId().equals(rec.getBaseId());
	    }
	    
	    public String toString() {
	        return Arrays.toString(lines);
	    }
	    
	    public int hashcode() {
	        return Arrays.hashCode(lines);
	    }
	    
	    public boolean equals(Object obj) {
	        FastqRecord that = (FastqRecord) obj;
	        return Arrays.equals(this.lines, that.lines);
	    }
	    
	    public void stripNonReadInfoInId() {
	        int idx = lines[0].indexOf(" ");
	        if (idx > 0) {
	            lines[0] = lines[0].substring(0, idx);
	        }
	    }
	    
	    public void appendToId(String suffix) {
	        if (!lines[0].endsWith(suffix)) {
	            lines[0] = lines[0] + suffix;
	        }
	    }
	    
	    String getQuality() {
	    	return lines[3];
	    }
	    
	    public String getSequence() {
	    	return lines[1];
	    }
	    
	    private void setQuality(String quality) {
	    	lines[3] = quality;
	    }
	}

	
	public static void main(String[] args) throws Exception {
		long s = System.currentTimeMillis();

		String input = args[0];		
//		String input = "/home/lmose/code/abra/small.reads";
//		String input = "/home/lmose/code/abra/1.read";
//		String input = "/home/lmose/code/abra/10.reads";
//		String input = "/home/lmose/code/abra/100.fastq";
		
		
		new Mapzilla().run(input);
		
//		Mapzilla mapz = new Mapzilla();
//		mapz.initPositions();
//		mapz.findRefPosition(1233695019);
		
		long e = System.currentTimeMillis();
		
		System.err.println("Elapsed: " + (e-s));
	}
}
