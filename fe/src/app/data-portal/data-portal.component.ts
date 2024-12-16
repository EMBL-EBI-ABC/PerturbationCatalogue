import {Component, OnInit} from '@angular/core';
import {MatCardModule} from "@angular/material/card";
import {MatTableModule} from "@angular/material/table";
import {MatPaginator, MatPaginatorModule, PageEvent} from "@angular/material/paginator";
import {CommonModule} from "@angular/common";
import {FormsModule} from "@angular/forms";
import {MatFormFieldModule} from "@angular/material/form-field";
import {MatInputModule} from "@angular/material/input";
import {MatListModule} from "@angular/material/list";
import {ElasticService} from "../services/elastic.service";

@Component({
  selector: 'app-data-portal',
  standalone: true,
  imports: [
    MatCardModule,
    MatTableModule,
    MatPaginator,
    MatPaginatorModule,
    CommonModule,
    FormsModule,
    MatFormFieldModule,
    MatInputModule,
    MatCardModule,
    MatListModule,
  ],
  templateUrl: './data-portal.component.html',
  styleUrl: './data-portal.component.css'
})
export class DataPortalComponent implements OnInit {
  data: any;
  aggregations: any;
  filters = {
    sequenceType: [],
    geneCategory: [],
    publicationYear: []
  }
  totalResults = 0;
  pageSize = 15;
  currentPage = 0;
  searchQuery: string = '';

  displayedColumns: string[] = ['urn', 'title', 'sequenceType', 'geneName', 'geneCategory', 'publicationYear',
    'numVariants'];

  constructor(private elasticService: ElasticService) {}

  ngOnInit() {
    this.fetchData();
  }

  getKeys(obj: any): string[] {
    return Object.keys(obj);
  }

  fetchData() {
    const start = this.currentPage * this.pageSize;
    const size = this.pageSize;
    const query = this.searchQuery;
    this.elasticService.getData(start, size, query, this.filters).subscribe(
      (response) => {
        this.data = response.results.map((row: Record<string, string>) => {
          const parsedRow: Record<string, string> = {};
          for (const key of this.getKeys(row)) {
            parsedRow[key] = row[key];
          }
          return parsedRow;
        });
        this.totalResults = response.total;
        this.aggregations = response.aggregations;
      },
      (error) => {
        console.error('Error fetching data:', error);
      }
    );
  }

  onPageChange(event: PageEvent) {
    this.pageSize = event.pageSize;
    this.currentPage = event.pageIndex;
    this.fetchData();
  }

  onSearchChange(): void {
    this.currentPage = 0;
    this.fetchData();
  }

  onFilterClick(filterKey: string, filterValue: string): void {
    // @ts-ignore
    const index = this.filters[filterKey].indexOf(filterValue);
    if (index > -1) {
      // @ts-ignore
      this.filters[filterKey].splice(index, 1)
    } else {
      // @ts-ignore
      this.filters[filterKey].push(filterValue);
    }
    this.fetchData();
  }
}
