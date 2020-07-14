#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "output_connector_data_transfer/output_connector_data.h"
#include "output_connector_data_transfer/output_connection.h"
#include "output_connector_data_transfer/global_connections_by_slot_name.h"
#include "output_writer/generic.h"
#include "utility/string_utility.h"

const int VARIABLES_LINE_LENGTH = 48;  // number of characters for the solver structure and variables, afterwards there will be connection lines

//! iterate over all nested solvers
void SolverStructureVisualizer::DiagramGenerator::
generateDiagramRecursion(std::stringstream &result, std::vector<std::vector<int>> &internalConnectionLinesCurrentSolver, std::vector<int> &slotLineNos,
                         int depth, bool isFirstChildSolver, bool isLastChildSolver)
{
  if (!currentSolver_)
    return;

  // prepare the output connections
  SolverStructureVisualizer::parseOutputConnection(currentSolver_);

  slotLineNos.clear();

  std::stringstream lineStart;
  for (int i = 0; i < depth; i++)
    lineStart << "│ ";

  result << lineStart.str();
  if (isLastChildSolver && currentSolver_->outputSlots.empty() && currentSolver_->outputConnections.empty())
    result << "└";
  else
    result << "├";
    //result << "┌";
  result << "── " << currentSolver_->name << "\n";

  if (currentSolver_->description != "")
  {
    result << lineStart.str() << "│   (\"" << currentSolver_->description << "\")" << "\n";
  }

  // print the output slots of the solver, store the line nos of each slot to slotLineNos
  if (!currentSolver_->outputSlots.empty())
  {
    result << lineStart.str() << "│  data slots: \n";

    // loop over the output connector slots of the current solver
    for (int i = 0; i < currentSolver_->outputSlots.size(); i++)
    {
      // print field variable and component name of the output slot
      std::string fieldVariableName = currentSolver_->outputSlots[i].fieldVariableName;
      std::string componentName = currentSolver_->outputSlots[i].componentName;
      std::string meshDescription = currentSolver_->outputSlots[i].meshDescription;
      std::string slotName = currentSolver_->outputSlots[i].slotName;

      std::stringstream s;
      s << lineStart.str() << "│  ";

      int meshReferenceNo = 0;

      // if mesh has already been references and thus is present in referencedMeshNames_
      if (std::find(referencedMeshNames_.begin(), referencedMeshNames_.end(), meshDescription) != referencedMeshNames_.end())
      {
        meshReferenceNo = std::find(referencedMeshNames_.begin(), referencedMeshNames_.end(), meshDescription) - referencedMeshNames_.begin();
      }
      else
      {
        referencedMeshNames_.push_back(meshDescription);
        meshReferenceNo = referencedMeshNames_.size() - 1;
      }
      s << "[" << std::string(1,char('a'+meshReferenceNo)) << "] ";

      // only add component name if it is different from the variable name
      if (fieldVariableName != componentName)
      {
        s << fieldVariableName;

        // only add component name if there are multiple components or, or the component is not not named "0" which is the standard
        if (currentSolver_->outputSlots[i].nComponents > 1 || componentName != "0")
          s << "." << componentName;
      }
      else
      {
        s << fieldVariableName;
      }

      // add note about variableNo
      if (currentSolver_->outputSlots[i].variableNo != 1)
        s << " (in variable" << currentSolver_->outputSlots[i].variableNo << ")";

      std::string outputSlotString = s.str();

      std::string mappingString;
      if (!currentSolver_->mappingsWithinSolver.empty())
      {
        for (const solver_t::OutputConnectionRepresentation &mapping : currentSolver_->mappingsWithinSolver)
        {
          int fromSlotNo = mapping.fromSlot;
          int toSlotNo = mapping.toSlot;

          // connection top to bottom
          if (fromSlotNo < toSlotNo)
          {
            if (i == fromSlotNo)
            {
              mappingString += "┌ ";
            }
            else if (i == toSlotNo)
            {
              mappingString += "└»";
            }
            else if (fromSlotNo < i && i < toSlotNo)
            {
              mappingString += "│ ";
            }
            else
            {
              mappingString += "  ";
            }
          }
          else
          {
            // connection bottom to top
            if (i == toSlotNo)
            {
              mappingString += "┌»";
            }
            else if (i == fromSlotNo)
            {
              mappingString += "└ ";
            }
            else if (toSlotNo < i && i < fromSlotNo)
            {
              mappingString += "│ ";
            }
            else
            {
              mappingString += "  ";
            }
          }

        }
      }
      int mappingStringLength = StringUtility::stringLength(mappingString);

      // shorten
      int currentLineLength = StringUtility::stringLength(outputSlotString)+2;
      const int requiredLineLength = VARIABLES_LINE_LENGTH - mappingStringLength;

      // shorten string if neccessary, add to stringstream s
      if (currentLineLength >= requiredLineLength)
      {
        outputSlotString = outputSlotString.substr(0, outputSlotString.length() - (currentLineLength - requiredLineLength));
        currentLineLength = StringUtility::stringLength(outputSlotString)+2;
        s.str("");
        s << outputSlotString;
      }

      // make the length of the slot name equal to 6 characters
      int slotNameInitialLength = StringUtility::stringLength(slotName);
      if (slotNameInitialLength > 6)
      {
        slotName = slotName.substr(0, 6);
      }
      else if (slotNameInitialLength < 6)
      {
        slotName = std::string(6-slotNameInitialLength, ' ') + slotName;
      }

      s << std::string(requiredLineLength - currentLineLength, ' ')     // fill with spaces
        << mappingString      // add mapping string, if any
        << " " << slotName
        << "¤" << i << "\n";
      result << s.str();

      // store line no
      std::string currentResult = result.str();
      int currentLineNo = std::count(currentResult.begin(), currentResult.end(), '\n')-1;

      //LOG(DEBUG) << "current Result [" << currentResult << "] contains " << currentLineNo << " line breaks";
      slotLineNos.push_back(currentLineNo);
      meshDescriptionsOfSlots_[currentLineNo] = meshDescription;
    }

    result << lineStart.str() << "│\n";
  }

  // print the output slot connections of the solver
  std::vector<solver_t::OutputConnectionRepresentation> &outputConnections = currentSolver_->outputConnections;
  if (!outputConnections.empty())
  {
    result << lineStart.str() << "│  slot connections: \n";
    for (int i = 0; i < outputConnections.size(); i++)
    {
      if (outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ba)
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].toSlot << "¤ <─ ¤" << outputConnections[i].fromSlot << "\n";
      }
      else
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].fromSlot << "¤";

        switch (outputConnections[i].type)
        {
        case solver_t::OutputConnectionRepresentation::output_connection_t::ba:
        case solver_t::OutputConnectionRepresentation::output_connection_t::ab:
          result << " ─> ";
          break;
        case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse:
          result << " <═> ";  // or use ⇔ ?
          break;
        case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy:
          result << " <─> ";
          break;

        };

        result << "¤" << outputConnections[i].toSlot << "\n";
      }
    }
    result << lineStart.str() << "│\n";
  }

  // get the line nos of the slots for nested solvers 0 and 1
  std::vector<std::shared_ptr<solver_t>> &children = currentSolver_->children;
  std::vector<std::vector<int>> slotLineNosTerm(children.size());
  std::vector<std::vector<std::vector<int>>> internalConnectionLinesCurrentSolverTerm(children.size());

  std::shared_ptr<solver_t> ownSolver = currentSolver_;

  int i = 0;
  // loop over child solvers
  for (std::vector<std::shared_ptr<solver_t>>::iterator iter = children.begin();
       iter != children.end(); iter++, i++)
  {
    currentSolver_ = (*iter);

    bool isFirstChildSolver = i == 0;
    bool isLastChildSolver = i == children.size()-1;

    // recursively call generateDiagram on subsolvers
    generateDiagramRecursion(result, internalConnectionLinesCurrentSolverTerm[i], slotLineNosTerm[i], depth+1, isFirstChildSolver, isLastChildSolver);

    VLOG(1) << std::string(depth*2, ' ') << "child " << i << ", i:" << internalConnectionLinesCurrentSolverTerm[i];

    // add end of subsolver
    if (isLastChildSolver)
      result << lineStart.str() << "└\n";
    else
      result << lineStart.str() << "│\n";

    // internal connection lines of the the 2nd and further solvers are finished, the one of the first are passed on to the parent solver
    if ((isFirstChildSolver && ownSolver->hasInternalConnectionToFirstNestedSolver)
        || (children.size() == 2 && isLastChildSolver && ownSolver->hasInternalConnectionToSecondNestedSolver))
    {
      VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") first child solver gives " << internalConnectionLinesCurrentSolverTerm[i] << " own are: " << internalConnectionLinesCurrentSolver;
      internalConnectionLinesCurrentSolver.insert(internalConnectionLinesCurrentSolver.end(), internalConnectionLinesCurrentSolverTerm[i].begin(), internalConnectionLinesCurrentSolverTerm[i].end());
      VLOG(1) << std::string(depth*2, ' ') << "now own: " << internalConnectionLinesCurrentSolver;
    }
    else 
    {
      VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") " << i << ". child solver gives " << internalConnectionLinesCurrentSolverTerm[i] << " all are: " << internalConnectionLinesAll_;
      internalConnectionLinesAll_.insert(internalConnectionLinesAll_.end(), internalConnectionLinesCurrentSolverTerm[i].begin(), internalConnectionLinesCurrentSolverTerm[i].end());
      VLOG(1) << std::string(depth*2, ' ') << "now all: " << internalConnectionLinesAll_;
    }
  }

  // loop over the output connector slots of the current solver
  for (int i = 0; i < ownSolver->outputSlots.size(); i++)
  {
    // add slot to internal connection lines
    if (i >= internalConnectionLinesCurrentSolver.size())
    {
      internalConnectionLinesCurrentSolver.resize(i+1);
    }
    internalConnectionLinesCurrentSolver[i].push_back(slotLineNos[i]);
    VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") add slot " << slotLineNos[i] << ", now i: " << internalConnectionLinesCurrentSolver;
  }

  VLOG(1) << std::string(depth*2, ' ') << "internalConnectionLinesCurrentSolver:" << internalConnectionLinesCurrentSolver;
  VLOG(1) << std::string(depth*2, ' ') << "internalConnectionLinesAll_:" << internalConnectionLinesAll_;

  // add connections between slots of nested solvers to externalConnectionLines
  for (int i = 0; i < outputConnections.size(); i++)
  {
    if (slotLineNosTerm.size() >= 2)
    {
      if ((outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ab
        || outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy
        || outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
        && slotLineNosTerm[0].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[1].size() > outputConnections[i].toSlot
      )
      {
        ExternalConnectionLine newExternalConnectionLine;
        newExternalConnectionLine.lineNoFrom = slotLineNosTerm[0][outputConnections[i].fromSlot];
        newExternalConnectionLine.lineNoTo   = slotLineNosTerm[1][outputConnections[i].toSlot];
        newExternalConnectionLine.lineColumn = 0;
        newExternalConnectionLine.lineType = outputConnections[i].type;

        // if the meshes of the connected slots are different, it is a mapping
        newExternalConnectionLine.involvesMapping = meshDescriptionsOfSlots_[newExternalConnectionLine.lineNoFrom] != meshDescriptionsOfSlots_[newExternalConnectionLine.lineNoTo];

        externalConnectionLines_.push_back(newExternalConnectionLine);
      }
      else if (outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ba
        && slotLineNosTerm[1].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[0].size() > outputConnections[i].toSlot
      )
      {
        ExternalConnectionLine newExternalConnectionLine;
        newExternalConnectionLine.lineNoFrom = slotLineNosTerm[0][outputConnections[i].toSlot];
        newExternalConnectionLine.lineNoTo   = slotLineNosTerm[1][outputConnections[i].fromSlot];
        newExternalConnectionLine.lineColumn = 0;
        newExternalConnectionLine.lineType = outputConnections[i].type;

        // if the meshes of the connected slots are different, it is a mapping
        newExternalConnectionLine.involvesMapping = meshDescriptionsOfSlots_[newExternalConnectionLine.lineNoFrom] != meshDescriptionsOfSlots_[newExternalConnectionLine.lineNoTo];

        externalConnectionLines_.push_back(newExternalConnectionLine);
      }
      else
      {
        VLOG(1) << "slotLineNosTerm is not good. fromSlot: " << outputConnections[i].fromSlot
          << ", toSlot: " << outputConnections[i].toSlot << ", slot line nos term0: " << slotLineNosTerm[0]
          << ", term1: " << slotLineNosTerm[1];
      }
    }
    else
    {
      VLOG(1) << "slotLineNosTerm is too short: " << slotLineNosTerm.size() << ", but there are " << outputConnections.size() << " outputConnections.";
    }
  }

  // if the own solver does not have output slots but it has only a single child, use their output slots,
  // this is true for MultipleInstances
  if (slotLineNos.empty() && children.size() == 1)
  {
    slotLineNos = slotLineNosTerm[0];
  }
  slotLineNosAll_.insert(slotLineNosAll_.end(), slotLineNos.begin(), slotLineNos.end());
}

//! initialize the diagram generator, afterwards, generateDiagram can be called
void SolverStructureVisualizer::DiagramGenerator::
initialize(std::shared_ptr<solver_t> solverRoot)
{
  currentSolver_ = solverRoot;
}

std::string SolverStructureVisualizer::DiagramGenerator::
generateDiagramWithoutConnectionLines()
{
  std::stringstream diagramWithoutConnectionLines;

  std::vector<std::vector<int>> internalConnectionLinesCurrentSolver;
  std::vector<int> slotLineNos;

  // call the implementation 
  generateDiagramRecursion(diagramWithoutConnectionLines, internalConnectionLinesCurrentSolver, slotLineNos, 0, false, false);

  // add the rest of the internal connection lines that are still storade under internalConnectionLinesCurrentSolverTerm to internalConnectionLinesAll
  internalConnectionLinesAll_.insert(internalConnectionLinesAll_.end(), internalConnectionLinesCurrentSolver.begin(), internalConnectionLinesCurrentSolver.end());

  return std::string(diagramWithoutConnectionLines.str());
}

std::string SolverStructureVisualizer::DiagramGenerator::
generateFinalDiagram()
{
  std::string diagram;

  //! externalConnectionLines_ contains the following information: <lineNoFrom, lineNoTo, lineColumn, lineType>
  // lineColumn is the horizontal position of the vertical data connection line

  // call generateDiagram on the nested solvers, this also fills the variables
  // externalConnectionLines_, internalConnectionLinesAll_ and slotLineNosAll_
  std::string diagramWithoutConnectionLinesString = generateDiagramWithoutConnectionLines();

  // set the "lineColumn" entry in externalConnectionLines
  // loop over externalConnectionLines
  for (int i = 0; i < externalConnectionLines_.size(); i++)
  {
    int startLineNo = externalConnectionLines_[i].lineNoFrom;
    int endLineNo   = externalConnectionLines_[i].lineNoTo;
    int lineColumn  = 0;

    // count number of starting lines inside this line
    for (int j = 0; j < externalConnectionLines_.size(); j++)
    {
      if (i == j)
        continue;
      int startLineNo2 = externalConnectionLines_[j].lineNoFrom;
      int endLineNo2   = externalConnectionLines_[j].lineNoTo;

      if (startLineNo < startLineNo2 && startLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo == startLineNo2 && endLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo < startLineNo2 && endLineNo2 == endLineNo)
        lineColumn++;
    }
    externalConnectionLines_[i].lineColumn = lineColumn;
  }

  // erase and remove internal connections that consist of one slot only
  internalConnectionLinesAll_.erase(
    std::remove_if(internalConnectionLinesAll_.begin(), internalConnectionLinesAll_.end(), 
    [](std::vector<int> connection){
      return connection.size() == 1;
    }), 
    internalConnectionLinesAll_.end()
  );

  //LOG(DEBUG) << "externalConnectionLines: " << externalConnectionLines_;
  LOG(DEBUG) << "internalConnectionLines: " << internalConnectionLinesAll_;

  // determine number of columns that are needed for internal connections
  int nColumnsForInternalConnectionLines = 0;
  int nLinesDiagram = std::count(diagramWithoutConnectionLinesString.begin(), diagramWithoutConnectionLinesString.end(), '\n');
  LOG(DEBUG) << "nLinesDiagram: " << nLinesDiagram;
    
  std::map<int,int> internalConnectionLineColumns;   //< key = index into internalConnectionLinesAll_, value = columnNo of that line
  for (int currentLineNo = 0; currentLineNo < nLinesDiagram; currentLineNo++)
  {
    // determine active connection lines 
    std::set<int> activeLinesIndices;
    int activeLineIndex = 0;

    // loop over all connection lines that are present in the current row
    for (const std::vector<int> &connectionLine : internalConnectionLinesAll_)
    {
      int lineNoFirst = *std::min_element(connectionLine.begin(), connectionLine.end());
      int lineNoLast = *std::max_element(connectionLine.begin(), connectionLine.end());
      
      // if the currently considered connectionLine is active, i.e. present in the current line "currentLineNo"
      if (lineNoFirst <= currentLineNo && currentLineNo <= lineNoLast)
      {
        // ensure that the connectionLine gets a column assigned if it does not have one already
        
        // if the current activeLine has no assigned column yet
        if (internalConnectionLineColumns.find(activeLineIndex) == internalConnectionLineColumns.end())
        {
          VLOG(1) << "find new column for internal connection line with index " << activeLineIndex << ": " << internalConnectionLinesAll_[activeLineIndex];

          // loop over possible columns until a free column is found that can be assign to the current activeLine
          for (int columnNo = 0; ; columnNo++)
          {
            // check if current column is still free
            bool columnIsFree = true;
            for (std::pair<int,int> columnNoConnectionLine : internalConnectionLineColumns)
            {
              if (columnNoConnectionLine.second == columnNo)
              {
                VLOG(1) << "column " << columnNo << " is occupied by line " << columnNoConnectionLine.first << ": " << internalConnectionLinesAll_[columnNoConnectionLine.first];
                columnIsFree = false;
                break;
              }
            }
            if (columnIsFree)
            {
              VLOG(1) << "found free column " << columnNo;

              // assign column to this connection line
              internalConnectionLineColumns[activeLineIndex] = columnNo;
              nColumnsForInternalConnectionLines = std::max(nColumnsForInternalConnectionLines, columnNo+1);
              break;
            }
          }
        }
      }

      activeLineIndex++;
    }
  }

  LOG(DEBUG) << "internalConnectionLineColumns: " << internalConnectionLineColumns;
  LOG(DEBUG) << "nColumnsForInternalConnectionLines: " << nColumnsForInternalConnectionLines;

  std::size_t pos = 0;
  diagram += DihuContext::globalConnectionsBySlotName()->getDescriptionForDiagram();

  diagram += "Solver structure: \n\n";

  // print actual diagram, into std::string diagram
  // loop over lines
  for (int currentLineNo = 0; pos < diagramWithoutConnectionLinesString.length(); currentLineNo++)
  {
    // determine current line by iterating over substring between "\n"s
    int nextPos = diagramWithoutConnectionLinesString.find("\n",pos+1);
    std::string line;

    if (pos == 0)
      pos = -1;

    if (nextPos != std::string::npos)
    {
      line = diagramWithoutConnectionLinesString.substr(pos+1, nextPos-pos-1);
    }
    else
    {
      line = diagramWithoutConnectionLinesString.substr(pos+1);
    }

    pos = nextPos;

    // here, line contains the current line that should be handled without "\n"

    // fill current line with spaces such that it has the correct length
    int lineLength = StringUtility::stringLength(line);

    // if the line is too short, add space until is has the correct length
    for (int i = 0; i < VARIABLES_LINE_LENGTH - lineLength; i++)
    {
      line += " ";
    }

#if 1
    // add internal connection lines

    std::string internalConnectionLinesColumns;

    // they are saved in internalConnectionLinesAll_

    bool connectionLineConnectsToTheRightInCurrentRow = false;

    // determine number of vertical lines in current row
    for (int currentColumnNo = 0; currentColumnNo < nColumnsForInternalConnectionLines; currentColumnNo++)
    {
      bool connectionLineFoundAtCurrentColumn = false;
      for (std::pair<int,int> internalConnectionLineColumn : internalConnectionLineColumns)
      {
        int columnNo = internalConnectionLineColumn.second;

        if (currentColumnNo == columnNo)
        {
          VLOG(1) << "currentColumnNo " << currentColumnNo << ", connection line " << internalConnectionLineColumn.first << "/" << internalConnectionLinesAll_.size();

          // vector of lines that are connected to the same connection line
          const std::vector<int> &connectionLine = internalConnectionLinesAll_[internalConnectionLineColumn.first];
          
          // begin and end row of connection line
          int lineNoFirst = *std::min_element(connectionLine.begin(), connectionLine.end());
          int lineNoLast = *std::max_element(connectionLine.begin(), connectionLine.end());
          
          VLOG(1) << "line " << currentLineNo << ", connection [" << lineNoFirst << "," << lineNoLast << "]";

          if (lineNoFirst <= currentLineNo && currentLineNo <= lineNoLast)
          {
            // check if it is connected to the current row
            if (std::find(connectionLine.begin(), connectionLine.end(), currentLineNo) != connectionLine.end())
            {
              internalConnectionLinesColumns += "├";
              connectionLineConnectsToTheRightInCurrentRow = true;
            }
            else
            {
              if (connectionLineConnectsToTheRightInCurrentRow)
              {
                internalConnectionLinesColumns += "÷";
              }
              else
              {
                internalConnectionLinesColumns += ":";
              }
            }

            connectionLineFoundAtCurrentColumn = true;
          }
          break;
        }
      }

      if (!connectionLineFoundAtCurrentColumn)
      {
        if (connectionLineConnectsToTheRightInCurrentRow)
        {
          internalConnectionLinesColumns += "─";
        }
        else
        {
          internalConnectionLinesColumns += " ";
        }
      }
    }   // loop over columns for internal connection lines

    if (connectionLineConnectsToTheRightInCurrentRow)
      internalConnectionLinesColumns += "──";
    else
      internalConnectionLinesColumns += "  ";

    std::string firstPart;
    int posFirstPart = 0;
    // add character by character
    for(int i = 0; i < line.size(); i++)
    {
      firstPart += line[i];
      if (StringUtility::stringLength(firstPart) == VARIABLES_LINE_LENGTH-2)
      {
        posFirstPart = i;
        break;
      }
    }
    line = firstPart + std::string(" ") + internalConnectionLinesColumns + line.substr(posFirstPart+1);

    //line += internalConnectionLinesColumns;
#endif

    // prepare for external connection lines
    // find out vertial data slot connection lines that go through this output line
    struct VerticalLineInCurrentRow
    {
      int lineColumn;         // the column in rhe current row where the vertical lines passes through
      bool lineStartsHere;    // if this row is the start of the line
      bool lineEndsHere;      // if this row is the end of the line
      bool hasMhere;          // if the current row contains the "m" that specifies a mapping
      solver_t::OutputConnectionRepresentation::output_connection_t type;  // type of the line
    };

    std::vector<VerticalLineInCurrentRow> verticalLinesInCurrentRow;    //! <lineColumn,thisLineStartsHere,thisLineEndsHere,type>, lineColumn
    bool lineStartsHere = false;
    bool lineEndsHere = false;
    bool unconnectedSlotHere = false;
    solver_t::OutputConnectionRepresentation::output_connection_t lineType = solver_t::OutputConnectionRepresentation::output_connection_t::ab;
    int lineColumn = 0;
#if 0
    // get maximum lineColumn
    int maximumLinePosition = 0;
    for (int i = 0; i < externalConnectionLines_.size(); i++)
    {
      maximumLinePosition = std::max((int)(externalConnectionLines_[i].lineColumn), maximumLinePosition);
    }
#endif
    // loop over all externalConnectionLines_ and see if they pass through the current line
    // collect all those rows into data structure verticalLinesInCurrentRow
    for (int i = 0; i < externalConnectionLines_.size(); i++)
    {
      // start and end row of current connection line
      int startLineNo = externalConnectionLines_[i].lineNoFrom;
      int endLineNo   = externalConnectionLines_[i].lineNoTo;

      // if the connection line starts or ends in the current row, there will be a horizontal line
      if (currentLineNo == startLineNo)
      {
        lineType = externalConnectionLines_[i].lineType;
        lineStartsHere = true;
        lineColumn = externalConnectionLines_[i].lineColumn;

        // add settings for current line pass to vector
        VerticalLineInCurrentRow verticalLineInCurrentRow;
        verticalLineInCurrentRow.lineColumn = lineColumn;
        verticalLineInCurrentRow.lineStartsHere = true;
        verticalLineInCurrentRow.lineEndsHere = false;
        verticalLineInCurrentRow.hasMhere = false;
        verticalLineInCurrentRow.type = lineType;
        verticalLinesInCurrentRow.push_back(verticalLineInCurrentRow);
      }
      else if (currentLineNo == endLineNo)
      {
        lineType = externalConnectionLines_[i].lineType;
        lineEndsHere = true;
        lineColumn = externalConnectionLines_[i].lineColumn;

        // add settings for current line pass to vector
        VerticalLineInCurrentRow verticalLineInCurrentRow;
        verticalLineInCurrentRow.lineColumn = lineColumn;
        verticalLineInCurrentRow.lineStartsHere = false;
        verticalLineInCurrentRow.lineEndsHere = true;
        verticalLineInCurrentRow.hasMhere = false;
        verticalLineInCurrentRow.type = lineType;
        verticalLinesInCurrentRow.push_back(verticalLineInCurrentRow);
      }
      else if (startLineNo < currentLineNo && currentLineNo < endLineNo)
      {
        // vertical line passes through current row
        solver_t::OutputConnectionRepresentation::output_connection_t lineType = externalConnectionLines_[i].lineType;
        int lineColumn = externalConnectionLines_[i].lineColumn;

        bool hasMhere = false;
        // if it involes a mapping and has the "m" here
        if (externalConnectionLines_[i].involvesMapping && (int)(0.5*(startLineNo+endLineNo)) == currentLineNo)
          hasMhere = true;

        // add settings for current line pass to vector
        VerticalLineInCurrentRow verticalLineInCurrentRow;
        verticalLineInCurrentRow.lineColumn = lineColumn;
        verticalLineInCurrentRow.lineStartsHere = false;
        verticalLineInCurrentRow.lineEndsHere = false;
        verticalLineInCurrentRow.hasMhere = hasMhere;
        verticalLineInCurrentRow.type = lineType;
        verticalLinesInCurrentRow.push_back(verticalLineInCurrentRow);
      }
    }

    // sort all found lines according to their column
    std::sort(verticalLinesInCurrentRow.begin(), verticalLinesInCurrentRow.end(), [](
      const VerticalLineInCurrentRow &a,
      const VerticalLineInCurrentRow &b)
    {
      return a.lineColumn < b.lineColumn;
    });

    //VLOG(1) << "line " << currentLineNo << ", currentLineNo " << currentLineNo << ", line=[" << line << "] has vertical lines: " << verticalLinesInCurrentRow;

    // fill current line with spaces such that it has the correct length
    const int requiredLineLength = VARIABLES_LINE_LENGTH + nColumnsForInternalConnectionLines + 5;  //< line length for structure and variables
    lineLength = StringUtility::stringLength(line);

    // determine fill character for horizontal line
    std::string fillCharacter = " ";
    if (lineStartsHere || lineEndsHere)
    {
      if (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
      {
        fillCharacter = "═";
      }
      else
      {
        fillCharacter = "─";
      }
    }
    else if (std::find(slotLineNosAll_.begin(), slotLineNosAll_.end(), currentLineNo) != slotLineNosAll_.end())
    {
      unconnectedSlotHere = true;
    }

    // if the line is yet too short, add fillCharacter until is has the correct length
    if (lineLength < requiredLineLength)
    {
      for (int i = 0; i < requiredLineLength - lineLength; i++)
      {
        // starting tip of arrow if the line starts or ends in the current row and has the according type
        if (i == 0 && lineStartsHere && (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::ba
          || lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy))
        {
          line += "<";
        }
        else if (i == 0 && lineEndsHere && (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::ab
          || lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy))
        {
          line += "<";
        }
        else
        {
          line += fillCharacter;
        }
      }
    }

    // add external connection lines
    // ---------
    // loop over the vertical lines in the current row
    int linesIndex = 0;
    if (!verticalLinesInCurrentRow.empty())
    {
      // loop over columns in the current row, current column is requiredLineLength + 2*i, i.e. it advances always by 2 characters
      for (int i = 0; i <= verticalLinesInCurrentRow.back().lineColumn; i++)
      {

        // determine events where line starts or ends
        bool lineCrosses = (lineStartsHere || lineEndsHere) && i < lineColumn;
        bool lineCurveStart = lineStartsHere && i == lineColumn;
        bool lineCurveEnd = lineEndsHere && i == lineColumn;

        // determine space character which is a horizontal line if the connection starts or ends here
        std::string space = " ";
        if (lineCrosses || ( (lineStartsHere || lineEndsHere) && i <= lineColumn))
        {
          if (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
          {
            space = "═";
          }
          else
          {
            space = "─";
          }
        }
        else if (unconnectedSlotHere && i == 0)
        {
          space = "x";
        }

        // if there is a vertical line passing at the current position
        if (i == verticalLinesInCurrentRow[linesIndex].lineColumn)
        {
          // depending on the type of the vertical line choose the right sign
          switch (verticalLinesInCurrentRow[linesIndex].type)
          {
          case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse:
            line += space;
            if (lineCrosses
                && !verticalLinesInCurrentRow[linesIndex].lineStartsHere
                && !verticalLinesInCurrentRow[linesIndex].lineEndsHere)
            {
              line += "╬";
            }
            else if (lineCurveStart)
            {
              line += "╗";
            }
            else if (lineCurveEnd)
            {
              line += "╝";
            }
            else if (verticalLinesInCurrentRow[linesIndex].lineStartsHere)
            {
              line += "╦";
            }
            else if (verticalLinesInCurrentRow[linesIndex].lineEndsHere)
            {
              line += "╩";
            }
            else if (verticalLinesInCurrentRow[linesIndex].hasMhere)
            {
              // if the meshes of the connected slots are different, it is a mapping, add an "m" in the center of the line
              line += "m";
            }
            else
            {
              line += "║";
            }
            break;
          case solver_t::OutputConnectionRepresentation::output_connection_t::ab:
          case solver_t::OutputConnectionRepresentation::output_connection_t::ba:
          case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy:
            line += space;
            if (lineCrosses
                && !verticalLinesInCurrentRow[linesIndex].lineStartsHere
                && !verticalLinesInCurrentRow[linesIndex].lineEndsHere)
            {
              line += "┼";
            }
            else if (lineCurveStart)
            {
              line += "┐";
            }
            else if (lineCurveEnd)
            {
              line += "┘";
            }
            else if (verticalLinesInCurrentRow[linesIndex].lineStartsHere)
            {
              line += "┬";
            }
            else if (verticalLinesInCurrentRow[linesIndex].lineEndsHere)
            {
              line += "┴";
            }
            else if (verticalLinesInCurrentRow[linesIndex].hasMhere)
            {
              // if the meshes of the connected slots are different, it is a mapping, add an "m" in the center of the line
              line += "m";
            }
            else
            {
              line += "│";
            }
            break;
          };
          linesIndex++;
        }
        else if (unconnectedSlotHere && i == 0)
        {
          line += std::string("x ");
        }
        else
        {
          // there is no vertical line passing through the current position, add two spaces
          line += space;
          line += space;
        }
      }
    }
    else if (unconnectedSlotHere)
    {
      line += std::string("x");
    }

    // add endline at the end of the current line
    line += "\n";

    // add the code line to the diagram
    diagram += line;
  }

  diagram += "Connection Types:\n";
  diagram += "  +··+   Internal connection, no copy\n";
  diagram += "  ════   Reuse variable, no copy\n";
  diagram += "  ───>   Copy data in direction of arrow\n";
  diagram += "  ─m──   Mapping between different meshes\n";
  diagram += "\n";
  diagram += "Referenced Meshes:\n";

  std::stringstream s;
  for (int i = 0; i < referencedMeshNames_.size(); i++)
  {
    s << "  [" << std::string(1, char('a'+i)) << "] " << referencedMeshNames_[i] << "\n";
  }
  diagram += s.str();

  //LOG(DEBUG) << "solver structure without connections:\n" << diagramWithoutConnectionLines.str() << "\n";
  LOG(DEBUG) << "solver structure with connections:\n" << diagram << "\n";

  return diagram;
}

//! produce the resulting file
void SolverStructureVisualizer::
writeDiagramFile(std::string filename)
{
  // only produce file on rank 0
  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    if (!filename.empty())
    {
      std::string diagram = getDiagram();

      // output diagram to a text file
      std::ofstream file;
      OutputWriter::Generic::openFile(file, filename);

      if (!file.is_open())
      {
        LOG(FATAL) << "Could not write to file \"" << filename << "\".";
      }

      // output diagram text
      file << diagram;
      file.close();

      LOG(INFO) << "File \"" << filename << "\" written.";
    }
  }

  LOG(DEBUG) << "SolverStructureVisualizer::clear data";

  // clear all values, this is required for unit tests where multiple different setups are created after each other
  solverRoot_ = nullptr;
  currentSolver_ = nullptr;
}
