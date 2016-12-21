package au.edu.qimr.utility;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.qcmg.common.log.QLogger;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import jdk.nashorn.internal.runtime.regexp.joni.constants.NodeType;

public class XmlCompare {

	    public boolean diff( String xml1, String xml2, List<String> diffs ) throws Exception
	    {
	        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	        dbf.setNamespaceAware(true);
	        dbf.setCoalescing(true);
	        dbf.setIgnoringElementContentWhitespace(true);
	        dbf.setIgnoringComments(true);
	        DocumentBuilder db = dbf.newDocumentBuilder();

	        Document doc1 = db.parse(new ByteArrayInputStream(xml1.getBytes()));
	        Document doc2 = db.parse(new ByteArrayInputStream(xml2.getBytes()));

	        doc1.normalizeDocument();
	        doc2.normalizeDocument();

	        return diff( doc1, doc2, diffs );
	    }

	    /**
	     * Diff 2 nodes and put the diffs in the list 
	     */
	    public boolean diff( Node node1, Node node2, List<String> diffs ) throws Exception
	    {
	        if( diffNodeExists( node1, node2, diffs ) ) return true;
	        	       
	        boolean nodeTypeDiff = diffNodeType(node1, node2, diffs );
	        boolean nodeValueDiff = nodeTypeDiff ? true : diffNodeValue(node1, node2, diffs );
	        	                 
	        //check attr and children if node value same 
	        if( !nodeValueDiff && node1.getNodeType() != Node.ATTRIBUTE_NODE){ 
	        	diffs.add( "node:" + node1.getNodeName() ) ;
		        diffAttributes( node1, node2, diffs );
		        diffChildNodes( node1, node2, diffs );		        
	        }
	        
	        return diffs.size() > 0;
	    }

	    /**
	     * Diff the nodes
	     */
	    public boolean diffChildNodes( Node node1, Node node2, List<String> diffs ) throws Exception
	    {
	        //Sort by Name
	        Map<String,Node> children1 = new LinkedHashMap<String,Node>();      
	        for( Node child1 = node1.getFirstChild(); child1 != null; child1 = child1.getNextSibling() ){
	        	if(child1.getNodeType() == Node.ELEMENT_NODE )
	        		children1.put( child1.getNodeName(), child1 );
	        }
	        
	        
	        //Sort by Name
	        Map<String,Node> children2 = new LinkedHashMap<String,Node>();      
	        for( Node child2 = node2.getFirstChild(); child2!= null; child2 = child2.getNextSibling() ) {
	        	if(child2.getNodeType() == Node.ELEMENT_NODE )
	        		children2.put( child2.getNodeName(), child2 );
	        }

	        //Diff all the children1
	        for( Node child1 : children1.values() ) {
	            Node child2 = children2.remove( child1.getNodeName() );
	            diff( child1, child2, diffs );
	        }

	        //Diff all the children2 left over
	        for( Node child2 : children2.values() ) {
	            Node child1 = children1.get( child2.getNodeName() );
	            diff( child1, child2, diffs );
	        }

	        return diffs.size() > 0;
	    }

	    /**
	     * Diff the nodes
	     */
	    public boolean diffAttributes( Node node1, Node node2, List<String> diffs ) throws Exception
	    {        
	        //Sort by Name
	        NamedNodeMap nodeMap1 = node1.getAttributes();
	        Map<String,Node> attributes1 = new LinkedHashMap<String,Node>();        
	        for( int index = 0; nodeMap1 != null && index < nodeMap1.getLength(); index++ )
	        {
	            attributes1.put( nodeMap1.item(index).getNodeName(), nodeMap1.item(index) );
	        }

	        //Sort by Name
	        NamedNodeMap nodeMap2 = node2.getAttributes();
	        Map<String,Node> attributes2 = new LinkedHashMap<String,Node>();        
	        for( int index = 0; nodeMap2 != null && index < nodeMap2.getLength(); index++ )
	        {
	            attributes2.put( nodeMap2.item(index).getNodeName(), nodeMap2.item(index) );

	        }

	        //Diff all the attributes1
	        for( Node attribute1 : attributes1.values() )
	        {
	            Node attribute2 = attributes2.remove( attribute1.getNodeName() );
	            diff( attribute1, attribute2, diffs );
	        }

	        //Diff all the attributes2 left over
	        for( Node attribute2 : attributes2.values() )
	        {
	            Node attribute1 = attributes1.get( attribute2.getNodeName() );
	            diff( attribute1, attribute2, diffs );
	        }

	        return diffs.size() > 0;
	    }
	    /**
	     * Check that the nodes exist
	     */
	    public boolean diffNodeExists( Node node1, Node node2, List<String> diffs ) throws Exception
	    {
	        if( node1 == null && node2 == null )
	        {
	            diffs.add( getPath(node2) + ": not exist in both XML" );
	            return true;
	        }

	        if( node1 == null && node2 != null )
	        {
	            diffs.add( getPath(node2) + ": not exists in second XML");
	            return true;
	        }

	        if( node1 != null && node2 == null )
	        {
	            diffs.add( getPath(node1) + ": not exists in first XML");
	            return true;
	        }

	        return false;
	    }

	    /**
	     * Diff the Node Type
	     */
	    public boolean diffNodeType( Node node1, Node node2, List<String> diffs ) throws Exception
	    {       
	        if( node1.getNodeType() != node2.getNodeType() ) 
	        {
	            diffs.add( getPath(node1) + ":type " + node1.getNodeType() + " != " + node2.getNodeType() );
	            return true;
	        }

	        return false;
	    }

	    /**
	     * Diff the Node Value
	     */
	    public boolean diffNodeValue( Node node1, Node node2, List<String> diffs ) throws Exception
	    {       
	        if( node1.getNodeValue() == null && node2.getNodeValue() == null ) {  return false;  }
	        
	    	if( node1.getNodeValue() == null || node2.getNodeValue() == null ) {
	    		diffs.add( getPath(node1) + " with null value!");
	            return true;
	        }
	        
	        if( !node1.getNodeValue().equals( node2.getNodeValue() ) ) {
	            diffs.add( getPath(node1) + " with differnt value!" );
	            return true;
	        }

	        return false;
	    }


	    /**
	     * Get the node path
	     */
	    public String getPath( Node node ){
	        StringBuilder path = new StringBuilder();
	        	        
	        if(node.getNodeType() == Node.ATTRIBUTE_NODE)
	        	path.insert( 0, "Atrribute::"+ node.getNodeName() );	        
	        else{
		        do {           
		            path.insert(0, node.getNodeName() );
		            path.insert( 0, "/" );
		        }while( ( node = node.getParentNode() ) != null );
	        
	        }
	        return path.toString();
	    }
	    
	    
	    public static Document createDocumentFromFile(File fXmlFile) throws Exception {
	    		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	    		DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	    		Document doc = dBuilder.parse(fXmlFile);

	    		//optional, but recommended
	    		//read this - http://stackoverflow.com/questions/13786607/normalization-in-dom-parsing-with-java-how-does-it-work
	    		doc.getDocumentElement().normalize();

//	    		System.out.println("Root element :" + doc.getDocumentElement().getNodeName());
	    		
	    		return doc; 	    	
	    } 
	    
	    
	    public static void main(final String[] args) throws Exception {
		//	QLogger logger =  options.getLogger(args);		
			try{    
				File flog = new File(args[3]);	
				
 				Document doc1 = createDocumentFromFile(new File(args[1]) );
 				Document doc2 = createDocumentFromFile(new File(args[2]) );
 				List<String> diffs = new ArrayList<String>();
 				XmlCompare compare =new XmlCompare();
 				compare.diff( doc1, doc2, diffs );
 				
 				for(String line : diffs) System.out.println(line);
 				
 				
			}catch(Exception e){ 
	        	System.err.println(e.toString());       	 
	        	System.exit(1);
	        }	
	    	
	    }

}
